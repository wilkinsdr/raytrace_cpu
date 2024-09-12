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
	char default_par_filename[] = "../par/disc_source_return_angdist.par";
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
	double dangle = par_file.get_parameter<double>("dangle", 0.1);

	const double r_isco = kerr_isco<double>(spin, +1);

	int Nangle = static_cast<int>(2*M_PI / dangle);

	double *escape_alpha_dist, *escape_beta_dist, *escape_norm_angle_dist, *escape_az_angle_dist;
	double *return_alpha_dist, *return_beta_dist, *return_norm_angle_dist, *return_az_angle_dist;
	double *lost_alpha_dist, *lost_beta_dist, *lost_norm_angle_dist, *lost_az_angle_dist;
	double *total_norm_angle_dist, *total_az_angle_dist;

	escape_alpha_dist = new double[Nangle];
	escape_beta_dist = new double[Nangle];
	escape_norm_angle_dist = new double[Nangle];
	escape_az_angle_dist = new double[Nangle];

	return_alpha_dist = new double[Nangle];
	return_beta_dist = new double[Nangle];
	return_norm_angle_dist = new double[Nangle];
	return_az_angle_dist = new double[Nangle];

	lost_alpha_dist = new double[Nangle];
	lost_beta_dist = new double[Nangle];
	lost_norm_angle_dist = new double[Nangle];
	lost_az_angle_dist = new double[Nangle];
	
	total_norm_angle_dist = new double[Nangle];
	total_az_angle_dist = new double[Nangle];
	
	for(int i=0; i<Nangle; i++)
	{
		escape_alpha_dist[i] = 0;
		escape_beta_dist[i] = 0;
		escape_norm_angle_dist[i] = 0;
		escape_az_angle_dist[i] = 0;
		return_alpha_dist[i] = 0;
		return_beta_dist[i] = 0;
		return_norm_angle_dist[i] = 0;
		return_az_angle_dist[i] = 0;
		lost_alpha_dist[i] = 0;
		lost_beta_dist[i] = 0;
		lost_norm_angle_dist[i] = 0;
		lost_az_angle_dist[i] = 0;
		total_norm_angle_dist[i] = 0;
		total_az_angle_dist[i] = 0;
	}

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

	//raytrace_source.redshift_start();
	raytrace_source.run_raytrace(1.1*r_esc, M_PI_2, show_progress);
	raytrace_source.range_phi();
    //raytrace_source.redshift(-1);

	raytrace_source.map_results(steps, t, r, theta, phi, redshift);

	for (int ray = 0; ray < raytrace_source.get_count(); ray++)
	{
		if (steps[ray] > 0)
		{
			if(raytrace_source.ray_cosalpha(ray) == 0) continue;
			//if(abs(redshift[ray] - 1) < 0.05) continue;

			double alpha = acos(raytrace_source.ray_cosalpha(ray));
			double beta = raytrace_source.ray_beta(ray);
			double norm_angle = acos(abs(sin(alpha)*sin(beta)));
			double az_angle = acos(-raytrace_source.ray_cosalpha(ray) / sqrt(1 - sin(alpha)*sin(alpha)*sin(beta)*sin(beta))); // minus sign here so the angle is wrt the inward radial
			if(sin(alpha)*cos(beta) < 0) az_angle *= -1;

			int alpha_bin = static_cast<int>((alpha + M_PI)/dangle);
			int beta_bin = static_cast<int>((beta + M_PI)/dangle);
			int norm_angle_bin = static_cast<int>((norm_angle + M_PI)/dangle);
			int az_angle_bin = static_cast<int>((az_angle + M_PI) / dangle);

			double ray_weight = (plane_iso) ? abs(sin(alpha)*sin(beta)) : 1;
			if(limb)
				ray_weight *= 1 + 2.06*(abs(sin(alpha)*sin(beta)));

			ray_count += (weight_norm) ? ray_weight : 1;
			
			total_norm_angle_dist[norm_angle_bin] += 1; // only 1 here to maintain normalisation, just fix drop-outs
			total_az_angle_dist[az_angle_bin] += 1;

			if (theta[ray] >= M_PI_2 && r[ray] >= r_isco && r[ray] < r_disc)
			{
				//if(abs(r[ray] - source_r) > 0.1*source_r || abs(phi[ray] - source_phi) > 0.1)
				//{
					return_count += ray_weight;
					return_alpha_dist[alpha_bin] += ray_weight;
					return_beta_dist[beta_bin] += ray_weight;
					return_norm_angle_dist[norm_angle_bin] += ray_weight;
					return_az_angle_dist[az_angle_bin] += ray_weight;

//					if(az_angle > 1.3 && az_angle < 1.7)
//					{
//						cout << "Naughty ray: " << alpha << " " << beta << " " << steps[ray] << " " << r[ray] << " " << theta[ray] << " " << phi[ray] << endl;
//					}
				//}

			}
			else if(r[ray] > r_esc)
			{
				escape_count += ray_weight;
				escape_alpha_dist[alpha_bin] += ray_weight;
				escape_beta_dist[beta_bin] += ray_weight;
				escape_norm_angle_dist[norm_angle_bin] += ray_weight;
				escape_az_angle_dist[az_angle_bin] += ray_weight;
			}
			else if(r[ray] < r_isco)
			{
				lost_count += ray_weight;
				lost_alpha_dist[alpha_bin] += ray_weight;
				lost_beta_dist[beta_bin] += ray_weight;
				lost_norm_angle_dist[norm_angle_bin] += ray_weight;
				lost_az_angle_dist[az_angle_bin] += ray_weight;
			}
            else
            {
				cout << "Naughty ray: " << alpha << " " << beta << " " << steps[ray] << " " << r[ray] << " " << theta[ray] << " " << phi[ray] << endl;
            }
		}

	}
	
	for(int i=0; i<Nangle; i++)
    {
	    double az_dist_sum = escape_az_angle_dist[i] +
                                return_az_angle_dist[i] +
                                    lost_az_angle_dist[i];

        escape_az_angle_dist[i] /= az_dist_sum;
        return_az_angle_dist[i] /= az_dist_sum;
        lost_az_angle_dist[i] /= az_dist_sum;
    }

	cout << endl << "Escape: " << escape_count / ray_count << endl;
	cout << "Return: " << return_count / ray_count << endl;
	cout << "Lost:   " << lost_count / ray_count << endl << endl;

	TextOutput outfile((char*)out_filename.c_str());
	for(int i=0; i<Nangle; i++)
	{
		double angle = -1*M_PI + i*dangle;
		outfile << angle << escape_alpha_dist[i] << escape_beta_dist[i] << escape_norm_angle_dist[i] << escape_az_angle_dist[i]
				<< return_alpha_dist[i] << return_beta_dist[i] << return_norm_angle_dist[i] << return_az_angle_dist[i]
				<< lost_alpha_dist[i] << lost_beta_dist[i] << lost_norm_angle_dist[i] << lost_az_angle_dist[i] << endl;
	}
	outfile.close();

	return 0;
}
