//
// Created by drw on 14/01/2020.
//
#include <iostream>
#include <string>
using namespace std;

#include "mapper/mapper_pointsource.h"
#include "include/par_file.h"
#include "include/par_args.h"

int main(int argc, char** argv)
{
	ParameterArgs par_args(argc, argv);

	// parameter configuration file
	char default_par_filename[] = "../par/pointsource_mapper.par";
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
	double V = par_file.get_parameter<double>("V");
	double spin = (par_args.key_exists("--spin")) ? par_args.get_parameter<double>("--spin")
	                                              :  par_file.get_parameter<double>("spin");
	double cosalpha0 = par_file.get_parameter<double>("cosalpha0", -0.9999);
	double cosalphamax = par_file.get_parameter<double>("cosalphamax", 0.9999);
	double dcosalpha = par_file.get_parameter<double>("dcosalpha");
	double beta0 = par_file.get_parameter<double>("beta0", -1*M_PI);
	double betamax = par_file.get_parameter<double>("betamax", M_PI);
	double dbeta = par_file.get_parameter<double>("dbeta");
	double r0 = par_file.get_parameter<double>("r0");
	double rmax = par_file.get_parameter<double>("rmax");
	int Nr = par_file.get_parameter<int>("Nr");
	bool logbin_r = par_file.get_parameter<bool>("logbin_r", false);
	double thetamax = par_file.get_parameter<double>("thetamax", M_PI_2);
	int Ntheta = par_file.get_parameter<int>("Ntheta");
	int Nphi = par_file.get_parameter<int>("Nphi");
	double source_velocity = par_file.get_parameter<double>("source_velocity");
	int source_vel_mode = par_file.get_parameter<int>("source_vel_mode", 0);
	int source_motion = par_file.get_parameter<int>("source_motion", 1);
	int show_progress = par_file.get_parameter<int>("show_progress", 1);

	Mapper_PointSource mapper(source, V, spin, dcosalpha, dbeta, r0, rmax, Nr, Ntheta, Nphi, logbin_r, thetamax);
	mapper.set_motion(source_velocity, source_motion, source_vel_mode);

	mapper.RedshiftStart();
	mapper.run_map(rmax, show_progress);
	mapper.save_hdf((char*)out_filename.c_str());

	return 0;
}
