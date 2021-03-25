//
// Created by Dan Wilkins on 3/17/21.
//

#include <iostream>
#include <cmath>
using namespace std;

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf.h>

#include "include/text_output.h"

// some GSL workspaces for integration and root finding
gsl_integration_workspace* integration_workspace;
gsl_root_fsolver* root_solver;

struct pcyg_params
{
    // wind parameters
    double w0;
    double beta;
    double gamma;
    double vinf;
    double k;
    double rout;
    double T;
    double turb;

    // the following are for passing dummy variables to GSL functions
    double r;
    double p;
    double z;
    double v;
};

inline double sigma(double r, pcyg_params* params)
{
    //
    // function containing the derivative of the velocity field
    //
    const double w = params->w0 + (1. - params->w0) * pow(1. - 1./r, params->beta);
    const double dwdr = params->beta * (1. - params->w0) * pow(1. - 1. / r, params->beta - 1.) / (r * r);

    return r * dwdr / w - 1.;
}

inline double tau0(double r, pcyg_params* params)
{
    //
    // radial component of the optical depth
    //
    const double w = params->w0 + (1. - params->w0) * pow(1. - 1./r, params->beta);
    const double dens = (0.01 * params->rout * params->rout) / (r*r*w);
//    return params->k * dens * r / (params->vinf * w);
//    return params->k * r / (params->vinf * w);

    const double dwdr = params->beta * (1. - params->w0) * pow(1. - 1. / r, params->beta - 1.) / (r * r);
    const double tau1 = params->T * (params->beta + 1.) * pow(1. - params->w0, -1 - params->beta) * pow(1. - w, params->gamma);
    return tau1 * r * dwdr / w;
}

double escape_integrand(double mu, void* p)
{
    //
    // integrand for evaluating escape/penetration probability of photons
    //
    pcyg_params *params = (pcyg_params*) p;
    const double t0 = tau0(params->r, params);
    const double s = sigma(params->r, params);

    return ((1. + s*mu*mu) / t0) * (1 - exp(-1 * t0/(1. + s*mu*mu)));
}

double escape_prob(double r, pcyg_params* params)
{
    //
    // angle-averaged escape probability for line photons emitted at radius r
    // for use in source function
    //
    double result, error;

    params->r = r;

    gsl_function escape_prob_func;
    escape_prob_func.function = &escape_integrand;
    escape_prob_func.params = params;

    gsl_integration_qags(&escape_prob_func, 0, 1, 0, 1e-7, 1000, integration_workspace, &result, &error);

    return result;
}

double continuum_penetration_prob(double r, pcyg_params* params)
{
    //
    // continuum photon penetration probability to radius r in wind
    // for use in source function
    //
    double result, error;

    params->r = r;

    gsl_function escape_prob_func;
    escape_prob_func.function = &escape_integrand;
    escape_prob_func.params = params;

    gsl_integration_qags(&escape_prob_func, sqrt(1. - (1./(r*r))), 1, 0, 1e-7, 1000, integration_workspace, &result, &error);

    return 0.5 * result;
}

double los_vel_z(double z, void* p)
{
    //
    // equation to solve position z along line of sight with impact parameter p
    // for a given value of the line-of-sight velocity (w*mu)
    // to be used with root finder
    //
    pcyg_params *params = (pcyg_params*) p;
    const double w = params->w0 + (1. - params->w0) * pow(1. - 1./sqrt((params->p)*(params->p) + z*z), params->beta);
    return w * z /sqrt((params->p)*(params->p) + z*z) - params->v;
}

double find_los_z(double wmu, double p, pcyg_params* params)
{
    //
    // find roots of los_vel_z
    // to solve position z along line of sight with impact parameter p
    // for a given value of the line-of-sight velocity (w*mu)
    //
    int status;
    int iter = 0, max_iter = 100;
    //double x_lo = (p>1) ? -1*params->rout : sqrt(1. - p*p), x_hi = params->rout;
    double x_lo = -1*params->rout, x_hi = (p>1) ? params->rout : -1*sqrt(1. - p*p);
    double z;
    params->v = wmu;
    params->p = p;

    gsl_function los_vel_z_func;
    los_vel_z_func.function = &los_vel_z;
    los_vel_z_func.params = params;

    // check that the equation has roots in the specified range (i.e. within the wind)
    // otherwise, don't even try
    if((los_vel_z(x_lo, params)*los_vel_z(x_hi, params)) > 0) return NAN;

    gsl_root_fsolver_set(root_solver, &los_vel_z_func, x_lo, x_hi);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate(root_solver);
        z = gsl_root_fsolver_root(root_solver);
        x_lo = gsl_root_fsolver_x_lower(root_solver);
        x_hi = gsl_root_fsolver_x_upper(root_solver);
        status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);

        if(status == GSL_SUCCESS) return z;
    }
    while(status == GSL_CONTINUE && iter < max_iter);

    // if we get this far, root finder didn't converge
    return NAN;
}

double source_func(double r, pcyg_params* params)
{
//    return 0;
    return continuum_penetration_prob(r, params)/escape_prob(r, params);
//    return 0.5*(1. - sqrt(1. - 1./(r*r)));
//      return 0;
}

double tau(double z, void* par)
{
    //
    // optical depth to photons leaving position (p,z) in cloud
    // in Sobolev approximation (steep velocity gradient), resonance line photons are
    // only re-emitted close to their emission point, according to velocity gradient
    //
    pcyg_params* params = (pcyg_params*) par;
    const double p = params->p;

    const double z0 = find_los_z(params->v, p, params);
    const double mu = z0/sqrt(p*p + z0*z0);
    const double r = sqrt(p*p + z0*z0);

    const double w = params->w0 + (1. - params->w0) * pow(1. - 1./r, params->beta);
    const double dwdr = params->beta * (1. - params->w0) * pow(1. - 1. / r, params->beta - 1.) / (r * r);
    const double s = sigma(sqrt(p*p + z0*z0), params);
    const double dwzdz = (1 + mu*mu*s) * w / r;

//    const double profile = (1./sqrt(M_PI)) * abs(gsl_sf_erf((mu*w - params->v)/params->turb) - gsl_sf_erf((-1 - params->v)/params->turb));
//    return profile * tau0(sqrt(p*p + z0*z0), params) / (1. + sigma(p*p + z0*z0, params)*mu*mu);

    if(z >= z0) // don't reabsorb yourself!
        //return tau0(sqrt(p*p + z0*z0), params) / (1. + sigma(p*p + z0*z0, params)*mu*mu);
        return tau0(sqrt(p*p + z0*z0), params) / (1. + sigma(p*p + z0*z0, params)*mu*mu);
    else
        return 0;
}

double tau_integrand(double z, void* par)
{
    //
    // optical depth to photons leaving position (p,z) in cloud
    // in Sobolev approximation (steep velocity gradient), resonance line photons are
    // only re-emitted close to their emission point, according to velocity gradient
    //
    pcyg_params* params = (pcyg_params*) par;
    const double p = params->p;
    const double mu = z/sqrt(p*p + z*z);
    const double r = sqrt(p*p + z*z);

    const double w = params->w0 + (1. - params->w0) * pow(1. - 1./sqrt((params->p)*(params->p) + z*z), params->beta);
    const double dwdr = params->beta * (1. - params->w0) * pow(1. - 1. / r, params->beta - 1.) / (r * r);
    const double s = sigma(sqrt(p*p + z*z), params);
    const double dwzdz = (1 + mu*mu*s) * w / r;

    const double delta_w = mu*w - params->v;
    const double phi = (1./sqrt(M_PI)) * (1./params->turb) * exp(-1 * (delta_w / params->turb)*(delta_w / params->turb));

    return tau0(r, params) * phi * dwzdz / (1. + s*mu*mu);
}

double integrate_tau(double wmu, double zstop, pcyg_params* params)
{
    //
    // continuum photon penetration probability to radius r in wind
    // for use in source function
    //
    double result, error;

    params->v = wmu;

    gsl_function tau_func;
    tau_func.function = &tau_integrand;
    tau_func.params = params;

    gsl_integration_qags(&tau_func, -1*params->rout, zstop, 0, 1e-7, 1000, integration_workspace, &result, &error);

    return result;
}

double flux_integrand(double p, void* par)
{
    pcyg_params *params = (pcyg_params*) par;

    const double z = find_los_z(params->v, p, params);
    //cout << "z = " << z << " for wmu = " << params->v << endl;

    params->p = p;

    const double r = sqrt(p*p + z*z);
    const double mu = z / sqrt(p*p + z*z);
    const double s = sigma(sqrt(p*p + z*z), params);

    const double this_tau = tau(z, params);;
    const double tau_star = tau(sqrt(1. - p*p), params);
    if(p < 1)
    {
        if(isnan(z)) return 2*p;
        return 2 * p * (source_func(r, params) * exp(-1 * this_tau) * tau0(r, params) / (1. + s*mu*mu) + exp(-1 * tau_star));
    }
    else
    {
        if(isnan(z)) return 0;
        return 2 * p * (source_func(r, params) * exp(-1 * this_tau) * tau0(r, params) / (1. + s*mu*mu));
    }
}

double integrate_flux(double wmu, pcyg_params* params)
{
    //
    // continuum photon penetration probability to radius r in wind
    // for use in source function
    //
    double result1, result2, error;

    params->v = wmu;

    gsl_function flux_func;
    flux_func.function = &flux_integrand;
    flux_func.params = params;

    gsl_integration_qags(&flux_func, 0, 1, 0, 1e-4, 10000, integration_workspace, &result1, &error);
    gsl_integration_qags(&flux_func, 1, params->rout, 0, 1e-4, 10000, integration_workspace, &result2, &error);

    return result1;
}

double integrate_flux2(double wmu, pcyg_params* params)
{
    double integral = 0;
    double dp = 0.1;

    params->v = wmu;

    for(double p=0; p<params->rout; p+=dp)
    {
        integral += flux_integrand(p, params) * dp;
    }
    return integral / 1.1;
}

int main()
{
    integration_workspace = gsl_integration_workspace_alloc(1000000);

    const gsl_root_fsolver_type* solver_type;
    solver_type = gsl_root_fsolver_brent;
    root_solver = gsl_root_fsolver_alloc(solver_type);

    pcyg_params params;
    params.beta = 1;
    params.gamma = 1;
    params.w0 = 0.01;
    params.vinf = 0.35;
    params.k = 0.0001;
    params.rout = 100;
    params.T = 1;
    params.turb = 0.001;

    TextOutput outfile("../dat/pcyg_sei.dat");

//    for(double p=0; p<10; p+=0.1)
//    {
//        params.v = 0.9;
//        outfile << p << flux_integrand(p, &params) << endl;
//    }

    for(double wmu=-1.5; wmu<1.5; wmu+=0.01)
    {
        outfile << wmu << integrate_flux2(wmu, &params) << endl;
    }
    outfile.close();

    gsl_integration_workspace_free(integration_workspace);
    gsl_root_fsolver_free(root_solver);

    return 0;
}