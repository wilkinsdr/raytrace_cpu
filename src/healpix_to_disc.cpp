//
// Created by drw on 14/01/2020.
//
#include <iostream>
#include <string>
using namespace std;

#include "raytracer/healpix_pointsource.h"
#include "include/par_file.h"
#include "include/par_args.h"
#include "include/text_output.h"

template <typename T>
T powerlaw3(T r, T q1, T rb1, T q2, T rb2, T q3)
{
    if(r < rb1)
        return pow(r, -1*q1);
    else if (r < rb2)
        return pow(rb1, q2-q1) * pow(r, -1*q2);
    else
        return pow(rb1, q2-q1) * pow(rb2, q3-q2) * pow(r, -1*q3);
}

int main(int argc, char** argv)
{
	ParameterArgs par_args(argc, argv);

	int *steps;
	double *t, *r, *theta, *phi, *redshift;

	// parameter configuration file
	char default_par_filename[] = "../par/healpix_to_disc.par";
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
	int order = par_file.get_parameter<int>("order");
    int motion = par_file.get_parameter<int>("motion", 0);
    int basis = par_file.get_parameter<int>("basis", 1);
    double rmax = par_file.get_parameter<double>("rmax", 1000);
    double r_disc = par_file.get_parameter<double>("r_disc", 500);
	int show_progress = par_file.get_parameter<int>("show_progress", 1);
    double q_in = par_file.get_parameter<double>("q_in", 7);
    double rb1 = par_file.get_parameter<double>("rb1", 4);
    double q_mid = par_file.get_parameter<double>("q_mid", 0);
    double rb2 = par_file.get_parameter<double>("rb2", 10);
    double q_out = par_file.get_parameter<double>("q_out", 3);

	double r_isco = kerr_isco<double>(spin, +1);

	HealpixPointSource<double> raytrace_source(source, V, spin, order, motion, basis);

    raytrace_source.redshift_start(true);
    raytrace_source.run_raytrace(rmax, M_PI_2, show_progress);
    raytrace_source.range_phi();
    raytrace_source.redshift(-1, true);

    raytrace_source.map_results(steps, t, r, theta, phi, redshift);

    TextOutput outfile((char*)out_filename.c_str());
    for(int pix=0; pix< raytrace_source.get_num_pix(); pix++)
    {
        double pix_r = r[5*pix+4];
        double pix_theta = theta[5*pix+4];
        double emis = powerlaw3(r[5*pix+4], q_in, rb1, q_mid, rb2, q_out) * pow(redshift[5*pix+4], -3);
        if(pix_r > r_isco && pix_theta > (M_PI_2 - 1E-2) && pix_r < r_disc)
            outfile << pix << r[5*pix+4] << phi[5*pix+4] << redshift[5*pix+4] << emis << endl;
        else
            outfile << pix << "nan" << "nan" << "nan" << 0 << endl;
    }
    outfile.close();

	return 0;
}
