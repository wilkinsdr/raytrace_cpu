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
	char default_par_filename[] = "../par/pointsource_to_disc.par";
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
    double rmax = par_file.get_parameter<double>("rmax", 1000);
    double r_disc = par_file.get_parameter<double>("r_disc", 500);
    double r_min = par_file.get_parameter<double>("r_min", -1);
    int Nr = par_file.get_parameter<int>("Nr", 100);
    bool logbin_r = par_file.get_parameter<bool>("logbin_r", true);
    double gamma = par_file.get_parameter<double>("gamma", 2);
	int show_progress = par_file.get_parameter<int>("show_progress", 1);

	const double r_isco = kerr_isco<double>(spin, +1);

	if(r_min == -1) r_min = r_isco;

	const double dr = (logbin_r) ? exp(log(r_disc / r_min) / (Nr)) : (r_disc - r_min) / (Nr);

	const int num_rays = ((cosalphamax - cosalpha0) / dcosalpha) * ((betamax - beta0) / dbeta);

	disc_r = new double[Nr];
	area = new double[Nr];
	disc_ray_count = new int[Nr];
	disc_redshift = new double[Nr];
	disc_emis = new double[Nr];

	for(int ir=0; ir<Nr; ir++)
    {
	    disc_r[ir] = (logbin_r) ? r_min * pow(dr, ir) : r_min + ir * dr;

	    const double annulus_dr = (logbin_r) ? disc_r[ir]*(dr - 1) : dr;
        area[ir] = integrate_disc_area(disc_r[ir], disc_r[ir] + annulus_dr, spin);
    }

	PointSource<double> raytrace_source(source, V, spin, TOL, dcosalpha, dbeta, cosalpha0, cosalphamax, beta0, betamax);

    EllipseDiscDestination<double>* my_destination = new EllipseDiscDestination<double>(r_disc, r_isco, major_axis, minor_axis);

    raytrace_source.redshift_start();
    //raytrace_source.run_raytrace(rmax, M_PI_2, show_progress);
    raytrace_source.run_raytrace(my_destination, 1.1 * dist); 
    raytrace_source.range_phi();
    //raytrace_source.redshift(-1);
    raytrace_source.redshift(my_destination, -1);

    raytrace_source.map_results(steps, t, r, theta, phi, redshift);

    for(int ray=0; ray<raytrace_source.get_count(); ray++)
    {
        //if(steps[ray] > 0 && theta[ray] > M_PI_2 - 1E-2 && r[ray] >= r_min && r[ray] < r_disc)
        if(steps[ray] > 0 && my_destination->stopping_fn(r, theta, phi, spin && r[ray] >= r_min && r[ray] < r_disc)
        {
            const int ir = (logbin_r) ? static_cast<int>( log(r[ray] / r_min) / log(dr)) : static_cast<int>((r[ray] - r_min) / dr);

            ++disc_ray_count[ir];
            disc_redshift[ir] += redshift[ray];
            disc_emis[ir] += 1./pow(redshift[ray], gamma);
        }
    }

    TextOutput outfile((char*)out_filename.c_str());
    for(int ir=0; ir<Nr; ir++)
    {
        disc_redshift[ir] /= disc_ray_count[ir];
        disc_emis[ir] /= area[ir];
        double ray_frac = static_cast<double>(disc_ray_count[ir]) / num_rays;
        double energy_frac = ray_frac / (area[ir] * pow(disc_redshift[ir], gamma));

        outfile << disc_r[ir] << area[ir] << ray_frac << disc_redshift[ir] << disc_emis[ir] << energy_frac << endl;
    }
    outfile.close();

	return 0;
}
