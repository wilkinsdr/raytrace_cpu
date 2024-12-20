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
#include "include/disc.h"


int main(int argc, char** argv)
{
	ParameterArgs par_args(argc, argv);

	int *steps;
	double *t, *r, *theta, *phi, *redshift;

	// parameter configuration file
	char default_par_filename[] = "../par/healpix_disc_source_photonfrac.par";
	char *par_filename;
	if (par_args.key_exists("--parfile"))
	{
		string par_filename_str = par_args.get_string_parameter("--parfile");
		par_filename = (char *) par_filename_str.c_str();
	} else
		par_filename = default_par_filename;

	double source[4];

	double *disc_r, *area, *disc_redshift, *disc_emis;
	int *disc_ray_count;

	ParameterFile par_file(par_filename);
//	string out_filename = (par_args.key_exists("--outfile")) ? par_args.get_parameter<string>("--outfile")
//	                                                         : par_file.get_parameter<string>("outfile");
	double source_r = (par_args.key_exists("--source_r")) ? par_args.get_parameter<double>("--source_r")
	                                              :  par_file.get_parameter<double>("source_r");
	double source_phi = par_file.get_parameter<double>("source_phi", 1.5707);
	double V = par_file.get_parameter<double>("V", -1);
	double spin = (par_args.key_exists("--spin")) ? par_args.get_parameter<double>("--spin")
	                                              :  par_file.get_parameter<double>("spin");
    int order = par_file.get_parameter<double>("order");
	int show_progress = par_file.get_parameter<int>("show_progress", 1);
	double r_esc = par_file.get_parameter<double>("r_esc", 1000);
	double r_disc = par_file.get_parameter<double>("r_esc", 500);

	const double r_isco = kerr_isco<double>(spin, +1);

	source[0] = 0.;
	source[1] = source_r;
	source[2] = M_PI_2 - 1E-6;
	source[3] = source_phi;

	if(V<0) V = disc_velocity<double>(source_r, spin, +1);

	double return_count = 0;
	double escape_count = 0;
	double lost_count = 0;
	double ray_count = 0;

	double this_beta0 = 0;
	double this_betamax = M_PI;

//	for(int beta_set=0; beta_set<2; beta_set++)
//	{
//		double this_beta0 = (beta_set == 0) ? -1*M_PI : M_PI/2;
//		double this_betamax = (beta_set == 0) ? -1*(M_PI/2) : M_PI;

		cout << "beta = " << this_beta0 << " to " << this_betamax << endl;

		HealpixPointSource<double> raytrace_source(source, V, spin, order, 0, 1);

		raytrace_source.set_disc_source();

//		raytrace_source.redshift_start();
		raytrace_source.run_raytrace(1.1*r_esc, M_PI_2, show_progress);
		raytrace_source.range_phi();
//		raytrace_source.redshift(-1);

		raytrace_source.map_results(steps, t, r, theta, phi, redshift);

		for(int pix=0; pix< raytrace_source.get_num_pix(); pix++)
		{
			int ray = 5*pix + 4;
			if (steps[ray] > 0)
			{
				++ray_count;

				if (theta[ray] >= M_PI_2 && r[ray] >= r_isco && r[ray] < r_disc)
				{
					if(abs(r[ray] - source_r) > 0.1*source_r || abs(phi[ray] - source_phi) > 0.1)
						++return_count;
				}
				else if(r[ray] > r_esc)
				{
					++escape_count;
				}
				else
				{
					++lost_count;
				}
			}

		}
	//}

	cout << endl << "Escape: " << escape_count / ray_count << endl;
	cout << "Return: " << return_count / ray_count << endl;
	cout << "Lost:   " << lost_count / ray_count << endl << endl;

	return 0;
}
