/*
 * trace_rays.cpp
 *
 *  Created on: 7 Mar 2019
 *      Author: drw
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

#include "sourcetracer/sourcetracer_imageplane.h"
#include "include/par_file.h"
#include "include/spectrum.h"
#include "include/array.h"


bool pointsource_stop(double t, double r, double theta, double phi, double* args)
{
	const double a = args[0];
	const double h = args[1];
	const double rad = args[2];

	const double x = sqrt(r*r + a*a)*sin(theta)*cos(phi);
	const double y = sqrt(r*r + a*a)*sin(theta)*sin(phi);
	const double z = r*cos(theta);

	if( (x*x + y*y + (z-h)*(z-h)) < rad ) return true;

	return false;
}

int main(int argc, char** argv)
{
	// parameter configuration file
	char default_par_filename[] = "../par/outflow_spectrum.par";
	char* par_filename;
	if(argc == 2)
		par_filename = argv[1];
	else
		par_filename = default_par_filename;


	SourceTracer_ImagePlane<double> *emis_source, *abs_source;

	ParameterFile par_file(par_filename);
	string out_filename = par_file.get_parameter<string>("outfile");
    string line_filename = par_file.get_parameter<string>("linefile");
	string emis_specfile = par_file.get_parameter<string>("emis_specfile");
    string abs_specfile = par_file.get_parameter<string>("abs_specfile");
	double dist = par_file.get_parameter<double>("dist");
	double incl = par_file.get_parameter<double>("incl");
	double plane_phi0 = par_file.get_parameter<double>("plane_phi0", 0);
	double spin = par_file.get_parameter<double>("spin");
	double x0 = par_file.get_parameter<double>("x0");
	double xmax = par_file.get_parameter<double>("xmax");
	int Nx = par_file.get_parameter<int>("Nx");
	double y0 = par_file.get_parameter<double>("y0");
	double ymax = par_file.get_parameter<double>("ymax");
	int Ny = par_file.get_parameter<int>("Ny");
	double en0 = par_file.get_parameter<double>("en0");
	double enmax = par_file.get_parameter<double>("enmax");
	int Nen = par_file.get_parameter<int>("Nen");
	bool logbin_en = par_file.get_parameter<bool>("logbin_en", false);
    double spec_en0 = par_file.get_parameter<double>("spec_en0");
    double spec_enmax = par_file.get_parameter<double>("spec_enmax");
    int spec_Nen = par_file.get_parameter<int>("spec_Nen");
    bool logbin_specen = par_file.get_parameter<bool>("logbin_specen", false);
	double source_size_xy = par_file.get_parameter<double>("source_size_xy");
	double source_size_z = par_file.get_parameter<double>("source_size_z");
	double source_vel = par_file.get_parameter<double>("source_vel");
	int source_motion = par_file.get_parameter<int>("source_motion", 1);
	double pointsource_h = par_file.get_parameter<double>("pointsource_h");
	double pointsource_r = par_file.get_parameter<double>("pointsource_r");
	double tau = par_file.get_parameter<double>("tau");
	double tol = par_file.get_parameter<double>("tol", TOL);
    int show_progress = par_file.get_parameter<double>("show_progress", 1);

	double dx = (xmax - x0) / (Nx - 1);
	double dy = (ymax - y0) / (Ny - 1);

    double den = (logbin_en) ? exp(log(enmax/en0)/(Nen-1)) : (enmax - en0)/(Nen - 1);
	double spec_den = (logbin_specen) ? exp(log(spec_enmax/spec_en0)/(spec_Nen-1)) : (spec_enmax - spec_en0)/(spec_Nen - 1);

	Spectrum<double> abs_spectrum((char*)abs_specfile.c_str());
    Spectrum<double> emis_spectrum((char*)emis_specfile.c_str());

	cout << "*****" << endl;
	cout << "Image plane at  d = " << dist << " , incl = " << incl << endl;
	cout << "Spin a = " << spin << endl;
	cout << "*****" << endl << endl;

	// -- Absorption - only count rays that reach the corona

    double stopping_args[] = {spin, pointsource_h, pointsource_r};

    abs_source = new SourceTracer_ImagePlane<double>(dist, incl, x0, xmax, dx, y0, ymax, dy, spin, en0, enmax, Nen, logbin_en, tol, plane_phi0);
    abs_source->set_source(source_size_xy, source_size_z, source_vel, source_motion);
    abs_source->redshift_start();
    abs_source->set_stopping_fn(pointsource_stop, stopping_args);
    abs_source->run_source_trace( 1.5*dist, M_PI_2, show_progress );
    int* ray_status = abs_source->get_status();
    double *ray_t, *ray_r, *ray_theta, *ray_phi, *ray_redshift;
    int *ray_steps;
    abs_source->map_results(ray_steps, ray_t, ray_r, ray_theta, ray_phi, ray_redshift);

    double *spec_tau;
    spec_tau = new double[spec_Nen];
    for(int ien=0; ien<spec_Nen; ien++)
        spec_tau[ien] = 0;

    double max_tau = 0;

    int pointsource_count = 0;
    for(int ray=0; ray< abs_source->get_count(); ray++)
    {
        if(ray_status[ray] == 2)
        {
            ++pointsource_count;
            for(int ien=0; ien<abs_source->Nen; ien++)
            {
                double enshift = (logbin_en) ? en0 * pow(den, ien) : en0 + ien*den;
                for(int bin=0; bin<emis_spectrum.bins; bin++)
                {
                    const double this_en = emis_spectrum.energy[bin] * enshift;
                    const int spec_ien = (logbin_specen) ? static_cast<int>( log(this_en / spec_en0) / log(spec_den)) : static_cast<int>((this_en - spec_en0) / spec_den);
                    spec_tau[spec_ien] += abs_spectrum.counts[bin] * abs_source->absorb[ray][ien];
                }
            }
        }
    }

    cout << pointsource_count << " rays reached pointsource" << endl;

    // average the optical depth over all rays that reach the continuum source
    // and normalise the optical depth by the maximum optical depth at any energy
    for(int ien=0; ien<spec_Nen; ien++)
    {
        spec_tau[ien] /= pointsource_count;
        if (spec_tau[ien] >= max_tau)
            max_tau = spec_tau[ien];
    }
    for(int ien=0; ien<spec_Nen; ien++)
    {
        spec_tau[ien] *= (tau/max_tau);
    }

    // -- the absorbed continuum --

    double *obs_continuum;
    obs_continuum = new double[spec_Nen];

    for(int ien=0; ien<spec_Nen; ien++)
    {
        obs_continuum[ien] = exp(-1 * spec_tau[ien]);
    }

    delete abs_source;

    cout << "got here" << endl;


    // -- Emission --

    cout << "Computing emission..." << endl;

    cout << "constructing source" << endl;
	emis_source = new SourceTracer_ImagePlane<double>(dist, incl, x0, xmax, dx, y0, ymax, dy, spin, en0, enmax, Nen, logbin_en, tol, plane_phi0);
	cout << "setting source" << endl;
	emis_source->set_source(source_size_xy, source_size_z, source_vel, source_motion);
    cout << "calculating redshift stat" << endl;
	emis_source->redshift_start();
    cout << "running raytrace" << endl;
	emis_source->run_source_trace( 1.5*dist, M_PI_2, show_progress );

    cout << "allocating arrays" << endl;

	double *emis;
	double *line_emis;
	emis = new double[spec_Nen];
	line_emis = new double[Nen];
	for(int ien=0; ien<spec_Nen; ien++)
		emis[ien] = 0;
    for(int ien=0; ien<Nen; ien++)
        line_emis[ien] = 0;

    cout << "working through rays" << endl;

	for(int ray=0; ray< emis_source->get_count(); ray++)
	{
		for(int ien=0; ien<Nen; ien++)
		{
		    line_emis[ien] += emis_source->emis[ray][ien];// * (tau/max_tau);
		    double enshift = (logbin_en) ? en0 * pow(den, ien) : en0 + ien*den;
		    for(int bin=0; bin<emis_spectrum.bins; bin++)
            {
                const double this_en = emis_spectrum.energy[bin] * enshift;
                const int spec_ien = (logbin_specen) ? static_cast<int>( log(this_en / spec_en0) / log(spec_den)) : static_cast<int>((this_en - spec_en0) / spec_den);
                emis[spec_ien] += emis_spectrum.counts[bin] * emis_source->emis[ray][ien] * (tau/max_tau);
            }
		}
	}

	//delete emis_source;

    // -- write out the results --

    TextOutput linefile((const char*)line_filename.c_str());
    for(int ien=0; ien<spec_Nen; ien++)
    {
        double en = (logbin_en) ? en0 * pow(den, ien) : en0 + den * ien;
        linefile << en << line_emis[ien] << endl;
    }
    linefile.close();

	TextOutput outfile((const char*)out_filename.c_str());
	for(int ien=0; ien<spec_Nen; ien++)
	{
	    double en = (logbin_specen) ? spec_en0 * pow(spec_den, ien) : spec_en0 + spec_den * ien;
		outfile << en << emis[ien] << obs_continuum[ien] << endl;
	}
    outfile.close();

	cout << "Done" << endl;

	return 0;
}
