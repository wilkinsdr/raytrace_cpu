//
// Created by drw on 14/01/2020.
//
#include <iostream>
#include <string>
using namespace std;

#include "raytracer/pointsource.h"
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
	char default_par_filename[] = "../par/disc_source_photonfrac.par";
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
    double cosalpha0 = par_file.get_parameter<double>("cosalpha0", -0.995);
    double cosalphamax = par_file.get_parameter<double>("cosalphamax", 0.995);
    double dcosalpha = par_file.get_parameter<double>("dcosalpha");
    double beta0 = par_file.get_parameter<double>("beta0", -1*M_PI);
    double betamax = par_file.get_parameter<double>("betamax", M_PI);
    double dbeta = par_file.get_parameter<double>("dbeta");
	int show_progress = par_file.get_parameter<int>("show_progress", 1);
	double r_esc = par_file.get_parameter<double>("r_esc", 1000);
	double r_disc = par_file.get_parameter<double>("r_esc", 500);
	bool plane_iso = par_file.get_parameter<bool>("plane_iso", true);
	bool limb = par_file.get_parameter<bool>("limb", false);
	bool weight_norm = par_file.get_parameter<bool>("weight_norm", true);

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

	cout << "beta = " << this_beta0 << " to " << this_betamax << endl;

	PointSource<double> raytrace_source(source, V, spin, TOL, dcosalpha, dbeta, cosalpha0, cosalphamax, this_beta0,
		                                    this_betamax);

	raytrace_source.redshift_start();
	raytrace_source.run_raytrace(1.1*r_esc, M_PI_2, show_progress);
	raytrace_source.range_phi();
	raytrace_source.redshift(-1);

	raytrace_source.map_results(steps, t, r, theta, phi, redshift);

	for (int ray = 0; ray < raytrace_source.get_count(); ray++)
	{
		if (steps[ray] > 0)
		{
			double alpha = acos(raytrace_source.ray_cosalpha(ray));
			double beta = raytrace_source.ray_beta(ray);

			double ray_weight = (plane_iso) ? abs(sin(alpha)*sin(beta)) : 1;
			if(limb)
				ray_weight *= 1 + 2.06*(abs(sin(alpha)*sin(beta)));

			ray_count += (weight_norm) ? ray_weight : 1;

			if (theta[ray] >= M_PI_2 && r[ray] >= r_isco && r[ray] < r_disc)
			{
				if(abs(r[ray] - source_r) > 0.1*source_r || abs(phi[ray] - source_phi) > 0.1)
					return_count += ray_weight;
			}
			else if(r[ray] > r_esc)
			{
				escape_count += ray_weight;
			}
			else if(r[ray] < r_isco)
			{
				lost_count += ray_weight;
			}
		}

	}

	cout << endl << "Escape: " << escape_count / ray_count << endl;
	cout << "Return: " << return_count / ray_count << endl;
	cout << "Lost:   " << lost_count / ray_count << endl << endl;

	return 0;
}
