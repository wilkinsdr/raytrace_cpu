#include <iostream>
#include <cmath>
#include <valarray>
using namespace std;

#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>

#include "include/text_output.h"

typedef double Real;
typedef valarray<double> RealArray;

namespace discwind_model
{
    // precision for integral over impact parameter
    const double precision = 10;
    const double dp_min = 0.1;

    struct discwind_params
    {
        // model parameters - in order for XSPEC
        Real line_en;       // rest frame energy of resonance line
        Real vinf;          // terminal velocity of wind
        Real tau_tot;       // total optical depth
        Real wind_angle;    // starting angle of wind from polar axis
        Real incl;          // inclination of line of sight to polar axis
        Real turb;          // turbulent/thermal velocity width of line (fraction of vinf)
        Real beta;          // power law index in velocity profile
        Real alpha1;        // power law index 1 in optical depth profile
        Real alpha2;        // power law index 2 in optical depth profile
        Real w0;            // velocity of wind at photosphere (fraction of vinf)
        Real rout;          // outer radius of wind (innder radius is 1)
        Real z;             // redshift
        bool continuum;     // include absorption of continuum (assumed to come from star at r=1)
        bool line_emis;     // include line emission from wind
    };

    inline Real sigma(Real r, discwind_params* params)
    {
        //
        // function containing the derivative of the velocity field
        //
        const Real w = params->w0 + (1. - params->w0) * pow(1. - 1./r, params->beta);
        const Real dwdr = params->beta * (1. - params->w0) * pow(1. - 1./r, params->beta - 1.) / (r*r);

        return r * dwdr / w - 1.;
    }

    Real tau_integrand(Real w, void* p)
    {
        //
        // provide the integrand to normalise the optical depth
        //
        discwind_params* params = (discwind_params*) p;
        return pow(w, params->alpha1) * pow(1. - w, params->alpha2);
    }

    Real integrate_tau(discwind_params* params)
    {
        // allocate a workspace to perform the integral using GSL
        gsl_integration_workspace* integration_workspace = gsl_integration_workspace_alloc(1000);
        // and set up a GSL function to integrate
        gsl_function tau_func;
        tau_func.function = &tau_integrand;
        tau_func.params = params;

        Real integral, error;
        gsl_integration_qags(&tau_func, 0, 1, 0, 1e-7, 1000, integration_workspace, &integral, &error);

        gsl_integration_workspace_free(integration_workspace);

        return integral;
    }

    inline Real tau0(double r, discwind_params* params)
    {
        //
        // radial component of the optical depth
        //

        // remember the old values of the optical depth parameters
        // so only run the normalisation integrals if we have to
        static Real alpha1_last(-1), alpha2_last(-1);
        static Real tau_norm(1);

        if(params->alpha1 != alpha1_last || params->alpha2 != alpha2_last)
        {
            // only integrate the optical depth again if its parameters have changed
            alpha1_last = params->alpha1;
            alpha2_last = params->alpha2;
            tau_norm = integrate_tau(params);
        }

        const Real w = params->w0 + (1. - params->w0) * pow(1. - 1./r, params->beta);
        const Real dwdr = params->beta * (1. - params->w0) * pow(1. - 1./r, params->beta - 1.) / (r*r);
        return params->tau_tot * pow(w, params->alpha1) * pow(1. - w, params->alpha2) * r * dwdr / (w * tau_norm);
    }

    inline Real source_func(Real r, discwind_params* params)
    {
        //
        // Source function for line emission in the wind
        //
        if(params->line_emis && r>1)
            return 0.5*(1. - sqrt(1. - 1./(r*r)));  // approximation from Castor 1970
        else
            return 0;
    }

    struct los_vel_z_params
    {
        // parameters for los_vel_z for GSL root solver
        discwind_params* model_params;
        Real p;
        Real v;
    };

    Real los_vel_z(double z, void* p)
    {
        //
        // Equation to solve position z along line of sight with impact parameter p
        // for a given value of the line-of-sight velocity (w*mu)
        // to be used with root finder
        //
        los_vel_z_params *params = (los_vel_z_params*) p;
        const double w = params->model_params->w0 + (1. - params->model_params->w0) * pow(1. - 1./sqrt((params->p)*(params->p) + z*z), params->model_params->beta);
        return w * z /sqrt((params->p)*(params->p) + z*z) - params->v;
    }

    Real find_los_z(double v, double p, discwind_params* params)
    {
        //
        // Find roots of los_vel_z to solve position z along line of sight with impact parameter p
        // for a given value of the line-of-sight velocity (v = w*mu)
        //
        int status;
        int iter = 0, max_iter = 100;

        Real x_lo = -1*params->rout, x_hi = (p>1) ? params->rout : -1*sqrt(1. - p*p);
        Real z;

        los_vel_z_params fparam;
        fparam.model_params = params;
        fparam.v = v;
        fparam.p = p;

        // set up a GSL root solver
        // if we don't care about thread safety, we could just do this once when model is first called
        const gsl_root_fsolver_type* solver_type;
        solver_type = gsl_root_fsolver_brent;
        gsl_root_fsolver* root_solver = gsl_root_fsolver_alloc(solver_type);
        // and a GSL function to solve
        gsl_function los_vel_z_func;
        los_vel_z_func.function = &los_vel_z;
        los_vel_z_func.params = &fparam;

        // check that the equation has roots in the specified range (i.e. within the wind)
        // otherwise, don't even try
        if((los_vel_z(x_lo, &fparam)*los_vel_z(x_hi, &fparam)) > 0) return NAN;

        gsl_root_fsolver_set(root_solver, &los_vel_z_func, x_lo, x_hi);

        do
        {
            iter++;
            status = gsl_root_fsolver_iterate(root_solver);
            z = gsl_root_fsolver_root(root_solver);
            x_lo = gsl_root_fsolver_x_lower(root_solver);
            x_hi = gsl_root_fsolver_x_upper(root_solver);
            status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);

            if(status == GSL_SUCCESS)
            {
                gsl_root_fsolver_free(root_solver);
                return z;
            }
        }
        while(status == GSL_CONTINUE && iter < max_iter);

        gsl_root_fsolver_free(root_solver);
        // if we get this far, root finder didn't converge
        return NAN;
    }

    Real tau(Real z, Real p, Real phi, Real v, discwind_params* params)
    {
        //
        // Calculate the optical depth along a line of sight, impact parameter p, up to position z (starting at z=-inf)
        // for emission at energy corresponding to velocity v
        //

        // find the integrand at the peak of the line profile if line profile is narrow
        const Real los_z = find_los_z(v, p, params);
        // if the velocity is outside the wind, set the location to the edge
        const Real z0 = isnan(los_z) ?
                        ((v < -0.5) ?
                         -1 * sqrt(params->rout*params->rout - p*p) :
                         ((p >= 1 && v > 0.5) ? sqrt(params->rout*params->rout - p*p) : -1 * sqrt(1 - p*p)))
                                     : los_z;

        const Real mu = z0 / sqrt(p*p + z0*z0);

        // really want to check if in wind here, otherwise

        // los velocities at the beginning and end of the integration
        const Real r_in = sqrt(p*p + z*z);
        const Real mu_in = z / sqrt(p*p + z*z);
        const Real w_in = params->w0 + (1. - params->w0) * pow(1. - 1. / r_in, params->beta);
        const Real w_out = params->w0 + (1. - params->w0) * pow(1. - 1. / params->rout, params->beta);
        const Real mu_out = -1 * sqrt(params->rout * params->rout - p * p) / params->rout;

        // integration over Gaussian line profile to include turbulence
        const Real profile = 0.5 * (gsl_sf_erf((w_in * mu_in - v) / params->turb) -
                                    gsl_sf_erf((w_out * mu_out - v) / params->turb));

        // angle from disc normal to this point
        const Real costheta = (p*sin(phi)*sin(params->incl) - z0*cos(params->incl)) / sqrt(p*p + z0*z0);
        // is this point in the wind? (we require theta > theta_wind and theta<PI/2 so it's above the disc)
        const Real in_wind = (costheta < params->wind_angle && costheta > 0) ? 1. : 0;

        //cout << costheta << " " << params->wind_angle << endl;

        return in_wind * profile * tau0(sqrt(p*p + z0*z0), params) / (1. + sigma(p*p + z0*z0, params) * mu*mu);
    }

    struct flux_integrand_params
    {
        // parameters for flux_integrand in a structure so we can eventually integrate using GSL
        discwind_params* model_params;
        Real v;
    };

    Real flux_integrand(Real p, Real phi, void* par)
    {
        //
        // Calculate the flux emitted along line of sight, impact parameter p, of emission at energy corresponding
        // to velocity v (passed inside par struct)
        //

        flux_integrand_params *params = (flux_integrand_params*) par;

        // find the integrand at the peak of the line profile if line profile is narrow
        const Real los_z = find_los_z(params->v, p, params->model_params);
        // if the velocity is outside the wind, set the location to the edge
        const Real z0 = isnan(los_z) ?
                        ((params->v < -0.5) ?
                         -1 * sqrt(params->model_params->rout*params->model_params->rout - p*p) :
                         ((p >= 1 && params->v > 0.5) ? sqrt(params->model_params->rout*params->model_params->rout - p*p) : -1.1 * sqrt(1 - p*p)))
                                     : los_z;

        const Real r = sqrt(p*p + z0*z0);

        // optical depth to the star (i.e. the continuum source)
        const Real tau_star = tau(-1 * sqrt(1. - p*p), p, phi, params->v, params->model_params);
        // optical depth to the point at which line emission at this velocity peaks (either going to the star or the far edge of the cloud, depending on line of sight)
        const Real this_tau = (p<1) ? tau_star : tau(params->model_params->rout, p, phi, params->v, params->model_params);

        // term describing line emission from the wind (per unit continuum intensity)
        const Real emission = source_func(r, params->model_params) * (1. - exp(-1 * this_tau));

        // angle from disc normal to this point - need to know if we're above the disc plane for the continuum emission
        const Real costheta = (p*sin(phi)*sin(params->model_params->incl) + sqrt(1. - p*p)*cos(params->model_params->incl));

        // whether we add the continuum term depends on whether this line of sight ends on the star
        if(p < 1 && params->model_params->continuum && costheta > 0)
        {
            return p * (emission + exp(-1 * tau_star));
        }
        else
        {
            return p * emission;
        }
    }

    Real integrate_flux(Real v, discwind_params* params)
    {
        //
        // Simple integration of flux at energy corresponding to velocity v
        // (with the adaptive step size, this seems to be faster than using GSL integration)
        //
        Real integral = 0;
        Real dp = 0.1;
        Real dphi = 0.1;

        flux_integrand_params fparam;
        fparam.model_params = params;
        fparam.v = v;

        Real continuum_norm = 0;
        for(Real p=0; p<params->rout; p+=dp)
        {
            // adaptive integration step
            dp = p / precision;
            if(dp < dp_min) dp = dp_min;

            for(Real phi=-1*M_PI; phi<M_PI; phi+=dphi)
            {
                integral += flux_integrand(p, phi, &fparam) * dp * dphi;

                // angle from disc normal to this point - need to know if we're above the disc plane for the continuum emission
                const bool above_disc = true; //p*sin(phi)*cos(M_PI_2 - params->incl) - z0*sin(M_PI_2 - params->incl)) > 0;

                if(p < 1 && above_disc) continuum_norm += p * dp * dphi;
            }
        }

        return integral / continuum_norm;
    }

    void do_discwind(const RealArray& energyArray, const RealArray& params, RealArray& fluxArray)
    {
        //
        // main model function using parameter array from XSPEC
        //
        discwind_params par;
        par.line_en = params[0];
        par.vinf = params[1];
        par.tau_tot = params[2];
        par.turb = params[3];
        par.beta = params[4];
        par.alpha1 = params[5];
        par.alpha2 = params[6];
        par.w0 = params[7];
        par.rout = params[8];
        par.z = params[9];
        par.continuum = params[10];
        par.line_emis = params[11];

        for(size_t ien=0; ien<energyArray.size()-1; ien++)
        {
            // velocity corresponding to this energy
            const Real this_en = (1 + par.z) * 0.5*(energyArray[ien] + energyArray[ien+1]);
            const Real v = (par.line_en - this_en) / (par.line_en * par.vinf);

            fluxArray[ien] = integrate_flux(v, &par);
        }
    }

    void do_discwind(const RealArray& energyArray, discwind_params& par, RealArray& fluxArray)
    {
        //
        // overload of main model function to be passed a parameter structure
        //
        for(size_t ien=0; ien<energyArray.size()-1; ien++)
        {
            // velocity corresponding to this energy
            const Real this_en = (1 + par.z) * 0.5*(energyArray[ien] + energyArray[ien+1]);
            const Real v = (par.line_en - this_en) / (par.line_en * par.vinf);

            fluxArray[ien] = integrate_flux(v, &par);
        }
    }
}

int main()
{
    discwind_model::discwind_params par;
    par.line_en = 6.4;
    par.vinf = 0.3;
    par.tau_tot = 1;
    par.wind_angle = cos(60 * M_PI/180.);
    par.incl = 60 * M_PI/180.;
    par.turb = 0.1;
    par.beta = 1;
    par.alpha1 = 6.4;
    par.alpha2 = 6.4;
    par.w0 = 0.01;
    par.rout = 100;
    par.z = 0;
    par.continuum = true;
    par.line_emis = true;

    const Real en0 = 1;
    const Real enmax = 10;
    const int Nen = 100;
    const Real den = (enmax - en0) / Nen;

    RealArray energyArray(Nen);
    RealArray fluxArray(Nen-1);
    for(size_t ien=0; ien<energyArray.size(); ien++)
    {
        energyArray[ien] = en0 + ien*den;
    }

    discwind_model::do_discwind(energyArray, par, fluxArray);

    TextOutput outfile("../dat/discwind.dat");
    for(size_t ien=0; ien<energyArray.size()-1; ien++)
    {
        outfile << 0.5*(energyArray[ien] + energyArray[ien+1]) << fluxArray[ien] << endl;
    }
    outfile.close();

    return 0;
}