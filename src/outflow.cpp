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

#include "raytracer/sourcetracer_imageplane.h"
#include "include/par_file.h"

bool pointsource_stop(double t, double r, double theta, double phi, double* args)
{
	const double a = args[0];
	const double h = args[1];
	const double rad = args[2];

	const double x = sqrt(r*r + a*a)*sin(theta)*cos(phi);
	const double y = sqrt(r*r + a*a)*sin(theta)*sin(phi);
	const double z = r*cos(theta);

	if( (x*x + y*y + (z-h)*(z-h)) < rad ) return true;

	//if(r < 2) return true;

	return false;
}

int main(int argc, char** argv)
{
	// parameter configuration file
	char default_par_filename[] = "../par/outflow.par";
	char* par_filename;
	if(argc == 2)
		par_filename = argv[1];
	else
		par_filename = default_par_filename;


	SourceTracer_ImagePlane<double> *RaytraceSource;

	ParameterFile par_file(par_filename);
	string out_filename = par_file.get_parameter<string>("outfile");
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
	double source_size_xy = par_file.get_parameter<double>("source_size_xy");
	double source_size_z = par_file.get_parameter<double>("source_size_z");
	double source_vel = par_file.get_parameter<double>("source_vel");
	int source_motion = par_file.get_parameter<int>("source_motion", 1);
	double pointsource_h = par_file.get_parameter<double>("pointsource_h");
	double pointsource_r = par_file.get_parameter<double>("pointsource_r");
	double tau = par_file.get_parameter<double>("tau");
	double tol = par_file.get_parameter<double>("tol", TOL);

	double dx = (xmax - x0) / (Nx - 1);
	double dy = (ymax - y0) / (Ny - 1);

	cout << "*****" << endl;
	cout << "Image plane at  d = " << dist << " , incl = " << incl << endl;
	cout << "Spin a = " << spin << endl;
	cout << "*****" << endl << endl;

	RaytraceSource = new SourceTracer_ImagePlane<double>(dist, incl, x0, xmax, dx, y0, ymax, dy, spin, en0, enmax, Nen, logbin_en, tol, plane_phi0);
	RaytraceSource->set_source(source_size_xy, source_size_z, source_vel, source_motion);
	RaytraceSource->RedshiftStart();
	RaytraceSource->run_source_trace( 1.5*dist );

	double *en, *emis, *abs;

	emis = new double[Nen];
	for(int ien=0; ien<Nen; ien++)
		emis[ien] = 0;

	for(int ray=0; ray< RaytraceSource->get_count(); ray++)
	{
		for(int ien=0; ien<Nen; ien++)
		{
			emis[ien] += RaytraceSource->emis[ray][ien];
		}
	}

	delete RaytraceSource;

	double stopping_args[] = {spin, pointsource_h, pointsource_r};

	RaytraceSource = new SourceTracer_ImagePlane<double>(dist, incl, x0, xmax, dx, y0, ymax, dy, spin, en0, enmax, Nen, logbin_en, tol, plane_phi0);
	RaytraceSource->set_source(source_size_xy, source_size_z, source_vel, source_motion);
	RaytraceSource->RedshiftStart();
	RaytraceSource->set_stopping_fn(pointsource_stop, stopping_args);
	RaytraceSource->run_source_trace( 1.5*dist );
	int* ray_status = RaytraceSource->get_status();
	double *ray_t, *ray_r, *ray_theta, *ray_phi, *ray_redshift;
	int *ray_steps;
    RaytraceSource->map_results(ray_steps, ray_t, ray_r, ray_theta, ray_phi, ray_redshift);

	double *obs_continuum, *total_absorb;
	obs_continuum = new double[RaytraceSource->Nen];
    total_absorb = new double[RaytraceSource->Nen];
	double max_tau = 0;

	int pointsource_count = 0;
	for(int ray=0; ray< RaytraceSource->get_count(); ray++)
	{
		if(ray_status[ray] == 2)
		{
			++pointsource_count;
			for(int ien=0; ien<RaytraceSource->Nen; ien++)
			{
				if(RaytraceSource->absorb[ray][ien] > max_tau)
					max_tau = RaytraceSource->absorb[ray][ien];
			}
		}
	}

	cout << "max tau = " << max_tau << endl;
	for(int ray=0; ray< RaytraceSource->get_count(); ray++)
	{
		if(ray_status[ray] == 2)
		{
			for(int ien=0; ien<RaytraceSource->Nen; ien++)
			{
				if(RaytraceSource->absorb[ray][ien] >= 0) obs_continuum[ien] += exp(-1 * RaytraceSource->absorb[ray][ien]*(tau/max_tau));
				total_absorb[ien] += RaytraceSource->absorb[ray][ien]*(tau/max_tau);
			}
		}
//		else
//        {
//            for(int ien = 0; ien < Nen; ien++)
//            {
//                emis[ien] += RaytraceSource->emis[ray][ien];
//            }
//        }
	}

	cout << pointsource_count << " rays reached pointsource" << endl;

	TextOutput outfile((const char*)out_filename.c_str());
	for(int ien=0; ien<Nen; ien++)
	{
		outfile << RaytraceSource->energy[ien] << emis[ien] << obs_continuum[ien] << total_absorb[ien] << endl;
	}
    outfile.close();

	cout << "Done" << endl;

	return 0;
}
