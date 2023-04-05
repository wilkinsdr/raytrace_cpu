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
	char default_par_filename[] = "../par/disc_source_photonfrac_r.par";
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
	string out_filename = (par_args.key_exists("--outfile")) ? par_args.get_parameter<string>("--outfile")
	                                                         : par_file.get_parameter<string>("outfile");
	double source_phi = par_file.get_parameter<double>("source_phi", 1.5707);
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
    double r_min = (par_args.key_exists("--rmin")) ? par_args.get_parameter<double>("--rmin")
                                                   :  par_file.get_parameter<double>("rmin", -1);
    int Nr = (par_args.key_exists("--Nr")) ? par_args.get_parameter<int>("--Nr")
                                           : par_file.get_parameter<int>("Nr");
    bool logbin_r = par_file.get_parameter<bool>("logbin_r", false);
	double r_disc = par_file.get_parameter<double>("r_esc", 500);
	bool plane_iso = par_file.get_parameter<bool>("plane_iso", true);
	bool limb = par_file.get_parameter<bool>("limb", false);
	bool weight_norm = par_file.get_parameter<bool>("weight_norm", true);

	const double r_isco = kerr_isco<double>(spin, +1);

	if(r_min == -1) r_min = r_isco;
    double dr = (logbin_r) ? exp(log(r_disc / r_min) / (Nr)) : (r_disc - r_min) / (Nr);

	source[0] = 0.;
	source[2] = M_PI_2 - 1E-6;
	source[3] = source_phi;

	TextOutput outfile((char*)out_filename.c_str());

    PointSource<double> *raytrace_source;

    for(int ir=0; ir<Nr; ir++)
    {
        double source_r = (logbin_r) ? r_min*pow(dr,ir) : r_min + ir*dr;

        cout << endl << ir+1 << '/' << Nr << ": r = " << source_r << endl;

        source[1] = source_r;

        double V = disc_velocity<double>(source_r, spin, +1);

        double return_count = 0;
        double escape_count = 0;
        double lost_count = 0;
        double ray_count = 0;

        raytrace_source = new PointSource<double> (source, V, spin, TOL, dcosalpha, dbeta, cosalpha0, cosalphamax, 0, M_PI);

        raytrace_source->redshift_start();
        raytrace_source->run_raytrace(1.1 * r_esc, M_PI_2, show_progress);
        raytrace_source->range_phi();
        raytrace_source->redshift(-1);

        raytrace_source->map_results(steps, t, r, theta, phi, redshift);

        for (int ray = 0; ray < raytrace_source->get_count(); ray++)
        {
            if (steps[ray] > 0)
            {
                double alpha = acos(raytrace_source->ray_cosalpha(ray));
                double beta = raytrace_source->ray_beta(ray);

                double ray_weight = (plane_iso) ? abs(sin(alpha) * sin(beta)) : 1;
                if (limb)
                    ray_weight *= 1 + 2.06 * (abs(sin(alpha) * sin(beta)));

                ray_count += (weight_norm) ? ray_weight : 1;

                if (theta[ray] >= M_PI_2 && r[ray] >= r_isco && r[ray] < r_disc)
                {
                    if (abs(r[ray] - source_r) > 0.1 * source_r || abs(phi[ray] - source_phi) > 0.1)
                        return_count += ray_weight;
                }
                else if (r[ray] > r_esc)
                {
                    escape_count += ray_weight;
                }
                else if (r[ray] < r_isco)
                {
                    lost_count += ray_weight;
                }
            }

        }

        double esc_frac = escape_count / ray_count;
        double return_frac = return_count / ray_count;
        double lost_frac = lost_count / ray_count;

        outfile << source_r << esc_frac << return_frac << lost_frac << endl;

        delete raytrace_source;
    }

    outfile.close();

	return 0;
}
