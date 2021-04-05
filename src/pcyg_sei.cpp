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
    bool continuum;
    bool line_emis;

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
//    const double dens = (0.01 * params->rout * params->rout) / (r*r*w);
//    //return params->k * dens * r / (params->vinf * w);
//    return params->k * r / (params->vinf * w);



    const double dwdr = params->beta * (1. - params->w0) * pow(1. - 1. / r, params->beta - 1.) / (r * r);
    const double tau1 = params->T * (params->beta + 1.) * pow(1. - params->w0, -1 - params->beta) * pow(1. - w, params->gamma);

    return params->T * pow(w, params->gamma) * pow(1. - w, params->gamma) * r * dwdr / w;

    //return tau1 * r * dwdr / w;
    return tau1;
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
    if(params->line_emis && r > 1)
        return 0.5*(1. - sqrt(1. - 1./(r*r)));
      //return continuum_penetration_prob(r, params)/escape_prob(r, params);
    else
        return 0;
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

    // find the integrand at the peak of the line profile if line profile is narrow
    // if the velocity is outside the wind, set the location to the edge
    const double los_z = find_los_z(params->v, p, params);
    const double z0 = isnan(los_z) ?
                      ((params->v < -0.5) ?
                            -1 * sqrt(params->rout*params->rout - p*p) :
                            ((p >= 1 && params->v > 0.5) ? sqrt(params->rout*params->rout - p*p) : -1 * sqrt(1 - p*p)))
                 : los_z;

    const double mu = z0 / sqrt(p*p + z0*z0);
    const double r = sqrt(p*p + z0*z0);
    const double s = sigma(sqrt(p*p + z0*z0), params);

    const double w = params->w0 + (1. - params->w0) * pow(1. - 1./sqrt((params->p)*(params->p) + z0*z0), params->beta);
    const double dwzdz = (1 + mu*mu*s) * w / r;

    if(params->turb > 0)
    {
        // los velocities at the beginning and end of the integration
        const double r_in = sqrt(p * p + z * z);
        const double mu_in = z / sqrt(p * p + z * z);
        const double w_in = params->w0 + (1. - params->w0) * pow(1. - 1. / r_in, params->beta);
        const double w_out = params->w0 + (1. - params->w0) * pow(1. - 1. / params->rout, params->beta);
        const double mu_out = -1 * sqrt(params->rout * params->rout - p * p) / params->rout;

        // integration over Gaussian line profile to include turbulence
        const double profile = 0.5 * (gsl_sf_erf((w_in * mu_in - params->v) / params->turb) -
                                         gsl_sf_erf((w_out * mu_out - params->v) / params->turb));
//        const double profile = 0.5 * (gsl_sf_erf((params->v - -1 * mu_out) / params->turb) -
//                                      gsl_sf_erf((params->v - w_in*mu_in) / params->turb));

//        cout << 0.5 * (gsl_sf_erf((w_in * mu_in - -0.01) / params->turb) -
//                       gsl_sf_erf((-1 * mu_out - -0.01) / params->turb))
//                       << " " <<
//                              0.5 * (gsl_sf_erf((w_in * mu_in - 0.01) / params->turb) -
//                                     gsl_sf_erf((-1 * mu_out - 0.01) / params->turb))
//                                     <<endl;

        //return tau0(sqrt(p * p + z0 * z0), params) / (1. + sigma(p * p + z0 * z0, params) * mu * mu);

 //        return sigma(p * p + z0 * z0, params)*mu*mu;
//       return profile * tau0(sqrt(p * p + z0 * z0), params) / (1. + sigma(p * p + z0 * z0, params) * mu * mu);
//        return 1./ (1. + sigma(p * p + z0 * z0, params) * mu * mu);
//        return tau0(sqrt(p * p + z0 * z0), params);

 //       return profile;
        //return z0;
        return profile * tau0(sqrt(p * p + z0 * z0), params) / (1. + sigma(p * p + z0 * z0, params) * mu * mu);
        //return tau0(sqrt(p * p + z0 * z0), params) / (1. + sigma(p * p + z0 * z0, params) * mu * mu);
    }
    else
    {
        // if no turbulence, then the profile is a delta function
        if(z > z0) // don't reabsorb yourself!
            //return tau0(sqrt(p*p + z0*z0), params) / (1. + sigma(p*p + z0*z0, params)*mu*mu);
            return tau0(sqrt(p * p + z0 * z0), params) / (1. + sigma(p * p + z0 * z0, params) * mu * mu);
        else
            return 0;
    }
}

double flux_integrand(double p, void* par)
{
    pcyg_params *params = (pcyg_params*) par;

    // find the integrand at the peak of the line profile if line profile is narrow
    // if the velocity is outside the wind, set the location to the edge
    const double los_z = find_los_z(params->v, p, params);
    const double z0 = isnan(los_z) ?
                      ((params->v < -0.5) ?
                       -1 * sqrt(params->rout*params->rout - p*p) :
                       ((p >= 1 && params->v > 0.5) ? sqrt(params->rout*params->rout - p*p) : -1.1 * sqrt(1 - p*p)))
                                   : los_z;

    const double mu = z0 / sqrt(p*p + z0*z0);

    params->p = p;

    const double r = sqrt(p*p + z0*z0);
    const double s = sigma(sqrt(p*p + z0*z0), params);

    const double w = params->w0 + (1. - params->w0) * pow(1. - 1./sqrt((params->p)*(params->p) + z0*z0), params->beta);
    const double dwzdz = (1 + mu*mu*s) * w / r;

    const double tau_star = tau(-1 * sqrt(1. - p*p), params);
    const double this_tau = (p<1) ? tau_star : tau(params->rout, params);

    if(params->turb > 0)
    {
        // los velocities at the beginning and end of the integration
        const double r_in = (p < 1) ? 1. : params->rout;
        const double mu_in = (p < 1) ? -1 * sqrt(1. - p*p) : sqrt(params->rout * params->rout - p * p) / params->rout;
        const double w_in = (p < 1) ? params->w0 : 1;
        const double mu_out = -1 * sqrt(params->rout * params->rout - p * p) / params->rout;

        // integration over Gaussian line profile to include turbulence
        const double profile = 0.5 * (gsl_sf_erf((w_in * mu_in - params->v) / params->turb) -
                                         gsl_sf_erf((mu_out - params->v) / params->turb));
//        const double profile = 0.5 * (gsl_sf_erf((params->v - -1 * mu_out) / params->turb) -
//                                         gsl_sf_erf((params->v - w_in*mu_in) / params->turb));



//        const double source = (1. - exp( -1 * profile * tau0(r, params) / (1. + s * mu * mu)))
//                              * source_func(r, params);

  //      return source_func(r, params) * (1. - exp(-1 * this_tau));

        const double source = source_func(r, params) * (1. - exp(-1 * this_tau));

  //      return dwzdz;

        //return (1. - exp( -1 * profile * tau0(r, params) / (1. + s * mu * mu)));

        //return profile * tau0(r, params) / (1. + s * mu * mu);

        if(p < 1 && params->continuum)
        {
            return 2 * p * (source + exp(-1 * tau_star));
        }
        else
        {
            return 2 * p * source;
        }
    }
    else
    {
        if(p < 1 && params->continuum)
        {
            if(isnan(z0) || abs(params->v) > 1) return 2 * p;
            return 2 * p * (source_func(r, params) * tau0(r, params) / (1. + s * mu * mu) +
                            exp(-1 * tau_star));
        } else
        {
            if(isnan(z0) || abs(params->v) > 1) return 0;
            return 2 * p * (source_func(r, params) * tau0(r, params) / (1. + s * mu * mu));
        }
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

//    gsl_integration_qags(&flux_func, 0, 1, 0, 1e-4, 10000, integration_workspace, &result1, &error);
//    gsl_integration_qags(&flux_func, 1, params->rout, 0, 1e-4, 10000, integration_workspace, &result2, &error);

    result2 = 0;
    gsl_integration_qags(&flux_func, 0, params->rout, 0, 1e-5, 1000, integration_workspace, &result1, &error);

    return result1 + result2;
}

double integrate_flux2(double wmu, pcyg_params* params)
{
    double integral = 0;
    double dp = 0.1;
    const double precision = 10;
    const double dp_min = 0.1;

    params->v = wmu;

    for(double p=0; p<params->rout; p+=dp)
    {
        dp = p / precision;
        if(dp < dp_min) dp = dp_min;
        integral += flux_integrand(p, params) * dp;
    }
    return integral / 1.1;
}

int main()
{
    integration_workspace = gsl_integration_workspace_alloc(1000);

    const gsl_root_fsolver_type* solver_type;
    solver_type = gsl_root_fsolver_brent;
    root_solver = gsl_root_fsolver_alloc(solver_type);

    pcyg_params params;
    params.beta = 1;
    params.gamma = 6.4;
    params.w0 = 0.01;
    params.vinf = 0.01;
    params.k = 0.0001;
    params.rout = 100;
    params.T = 1500;
    params.turb = 0.1;
    params.continuum = true;
    params.line_emis = true;

    TextOutput outfile("../dat/pcyg_sei.dat");

//    for(double wmu = -1.5; wmu<=1.5; wmu += 0.1)
//        for(double p = 0; p<10; p+=0.1)
//        {
//            outfile << find_los_z(wmu, p, &params) << p << endl;
//        }

//    for(double wmu=-1.5; wmu<1.5; wmu+=0.01)
//    {
//        params.v = wmu;
//        params.p = 0.5;
//        outfile << wmu << tau(100, &params) << endl;
//    }
//
//    for(double wmu=-1.5; wmu<1.5; wmu+=0.01)
//    {
//        params.v = wmu;
//        const double p = 1.01;
//        outfile << wmu << flux_integrand(p, &params) << endl;
//    }

//    for(double r=1; r<100; r+=0.1)
//    {
//        outfile << r << source_func(r, &params) << endl;
//    }
//
    for(double wmu=-1.5; wmu<1.5; wmu+=0.05)
    {
        outfile << wmu << integrate_flux(wmu, &params) << endl;
    }

//    const double en0 = 7.4;
//
//    for(double en = 3; en < 15; en+=0.05)
//    {
//        const double wmu = (en0 - en) / (en0 * params.vinf);
//        outfile << en << integrate_flux2(wmu, &params) << endl;
//    }
    outfile.close();

    gsl_integration_workspace_free(integration_workspace);
    gsl_root_fsolver_free(root_solver);

    return 0;
}