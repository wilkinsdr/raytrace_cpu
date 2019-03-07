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

#include "raytracer/imageplane.h"
#include "include/par_file.h"

int main(int argc, char** argv)
{
	// parameter configuration file
	char default_par_filename[] = "../par/trace_rays_imageplane.par";
	char* par_filename;
	if(argc == 2)
		par_filename = argv[1];
	else
		par_filename = default_par_filename;


	ImagePlane<double> *RaytraceSource;

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
	double tol = par_file.get_parameter<double>("tol", TOL);
	double write_step = par_file.get_parameter<double>("write_step", 10);
	double write_rmin = par_file.get_parameter<double>("write_step", -1);
	double write_rmax = par_file.get_parameter<double>("write_step", -1);
	double theta_max = par_file.get_parameter<double>("thetamax", 0);

	double dx = (xmax - x0) / (Nx - 1);
	double dy = (ymax - y0) / (Ny - 1);

	cout << "*****" << endl;
	cout << "Image plane at  d = " << dist << " , incl = " << incl << endl;
	cout << "Spin a = " << spin << endl;
	cout << "*****" << endl << endl;

	TextOutput outfile((const char*)out_filename.c_str());

	RaytraceSource = new ImagePlane<double>(dist, incl, x0, xmax, dx, y0, ymax, dy, spin, tol, plane_phi0);

	RaytraceSource->RunRaytrace( 1.5*dist, theta_max, &outfile, write_step, write_rmax, write_rmin );

	delete RaytraceSource;

	outfile.Close();

	cout << "Done" << endl;

	return 0;
}

