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
	double t0 = par_file.get_parameter<double>("t0");
	double tmax = par_file.get_parameter<double>("tmax");
	int Nt = par_file.get_parameter<int>("Nt");
	double source_size_xy = par_file.get_parameter<double>("source_size_xy");
	double source_size_z = par_file.get_parameter<double>("source_size_z");
	double source_vel = par_file.get_parameter<double>("source_vel");
	int source_motion = par_file.get_parameter<int>("source_motion", 1);
	double tol = par_file.get_parameter<double>("tol", TOL);

	double dx = (xmax - x0) / (Nx - 1);
	double dy = (ymax - y0) / (Ny - 1);

	cout << "*****" << endl;
	cout << "Image plane at  d = " << dist << " , incl = " << incl << endl;
	cout << "Spin a = " << spin << endl;
	cout << "*****" << endl << endl;

	RaytraceSource = new SourceTracer_ImagePlane<double>(dist, incl, x0, xmax, dx, y0, ymax, dy, spin, en0, enmax, Nen, logbin_en, t0, tmax, Nt, tol, plane_phi0);
	RaytraceSource->set_source(source_size_xy, source_size_z, source_vel, source_motion);
	RaytraceSource->RedshiftStart();
	RaytraceSource->run_source_trace( 1.5*dist );

	double *en, *emis, *abs;

	emis = new double[Nen];
	for(int ien=0; ien<Nen; ien++)
		emis[ien] = 0;

	for(int ray=0; ray<RaytraceSource->GetCount(); ray++)
	{
		for(int ien=0; ien<Nen; ien++)
		{
			emis[ien] += RaytraceSource->emis[ray][ien];
		}
	}

	delete RaytraceSource;

	TextOutput outfile((const char*)out_filename.c_str());
	for(int ien=0; ien<Nen; ien++)
	{
		outfile << RaytraceSource->energy[ien] << emis[ien] << endl;
	}
	outfile.Close();

	cout << "Done" << endl;

	return 0;
}
