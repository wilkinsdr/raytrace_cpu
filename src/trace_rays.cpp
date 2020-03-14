/*
 * trace_rays.cpp
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

#include "raytracer/pointsource.h"
#include "include/par_file.h"
#include "include/par_args.h"

void read_par(const char* filename, char* outfile, double* source, double& V, double& a, double& dcosalpha, double& dbeta, double& cosalpha0, double& cosalphamax, double& beta0, double& betamax, double& r_max, double& theta_max, int& write_step);

int main(int argc, char** argv)
{
	ParameterArgs par_args(argc, argv);

	// parameter configuration file
	char default_par_filename[] = "../par/trace_rays.par";
	char *par_filename;
	if (par_args.key_exists("--parfile"))
	{
		string par_filename_str = par_args.get_string_parameter("--parfile");
		par_filename = (char *) par_filename_str.c_str();
	} else
		par_filename = default_par_filename;

	double source[4];


	ParameterFile par_file(par_filename);
	string out_filename = (par_args.key_exists("--outfile")) ? par_args.get_parameter<string>("--outfile")
	                                                         : par_file.get_parameter<string>("outfile");
	par_file.get_parameter_array("source", source, 4);
	double V = par_file.get_parameter<double>("V", -1);
	double spin = (par_args.key_exists("--spin")) ? par_args.get_parameter<double>("--spin")
	                                              :  par_file.get_parameter<double>("spin");
	double cosalpha0 = par_file.get_parameter<double>("cosalpha0", -0.995);
	double cosalphamax = par_file.get_parameter<double>("cosalphamax", 0.995);
	double dcosalpha = par_file.get_parameter<double>("dcosalpha");
	double beta0 = par_file.get_parameter<double>("beta0", -1*M_PI);
	double betamax = par_file.get_parameter<double>("betamax", M_PI);
	double dbeta = par_file.get_parameter<double>("dbeta");
	double r_max = par_file.get_parameter<double>("r_max", 100);
	double theta_max = par_file.get_parameter<double>("theta_max", M_PI_2);
	int show_progress = par_file.get_parameter<int>("show_progress", 1);
	double write_step = par_file.get_parameter<double>("write_step", 10);
	double write_rmin = par_file.get_parameter<double>("write_rmin", -1);
	double write_rmax = par_file.get_parameter<double>("write_rmax", -1);
	bool write_cartesian = par_file.get_parameter<double>("write_cartesian", true);

	if(V < 0) V = disc_velocity(source[1], spin, +1);

	PointSource<double> *RaytraceSource;

	cout << "*****" << endl;
	cout << "Source r = [" << source[0] << " , " << source[1] << " , " << source[2] << " , " << source[3] << "] mu" << endl;
	cout << "Source angular velocity V = " << V << " mu*c" << endl;
	cout << "Spin a = " << spin << endl;
	cout << "*****" << endl << endl;

	TextOutput outfile(out_filename);

	RaytraceSource = new PointSource<double>( source, V, spin, TOL, dcosalpha, dbeta, cosalpha0, cosalphamax, beta0, betamax );
	RaytraceSource->run_raytrace(r_max, theta_max, show_progress, &outfile, write_step, write_rmax, write_rmin, write_cartesian);

	delete RaytraceSource;

	outfile.close();

	cout << "Done" << endl;

	return 0;
}


