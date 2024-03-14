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

	// parameter configuration file
	char default_par_filename[] = "../par/emissivity.par";
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
	double V = par_file.get_parameter<double>("V", 0);
	double spin = (par_args.key_exists("--spin")) ? par_args.get_parameter<double>("--spin")
	                                              :  par_file.get_parameter<double>("spin");
    double cosalpha0 = par_file.get_parameter<double>("cosalpha0", -0.995);
    double cosalphamax = par_file.get_parameter<double>("cosalphamax", 0.995);
    double dcosalpha = par_file.get_parameter<double>("dcosalpha");
    double beta0 = par_file.get_parameter<double>("beta0", -1*M_PI);
    double betamax = par_file.get_parameter<double>("betamax", M_PI);
    double dbeta = par_file.get_parameter<double>("dbeta");
	int show_f = (par_args.key_exists("--show_progress")) ? par_args.get_parameter<int>("--show_progress")
	                                        : par_file.get_parameter<int>("show_progress", 1);
	double r_max = par_file.get_parameter<double>("r_esc", 1000);
    double r_min = (par_args.key_exists("--rmin")) ? par_args.get_parameter<double>("--rmin")
                                                   :  par_file.get_parameter<double>("rmin", -1);
    int Nr = (par_args.key_exists("--Nr")) ? par_args.get_parameter<int>("--Nr")
                                           : par_file.get_parameter<int>("Nr", 100);
    double r_disc = par_file.get_parameter<double>("r_esc", 500);
    bool logbin_r = par_file.get_parameter<bool>("logbin_r", true);
    double gamma = par_file.get_parameter<double>("gamma", 2);
    double r_angle_disc_dis = par_file.get_parameter<double>("r_angle_disc_dis");
    double dist = par_file.get_parameter<double>("dist");

    double efficiency = par_file.get_parameter<double>("efficiency");
    double edd = par_file.get_parameter<double>("edd"); // accretion rate

    double major_axis = par_file.get_parameter<double>("major_axis");
    double minor_axis = par_file.get_parameter<double>("minor_axis");
    double thetalim = par_file.get_parameter<double>("thetalim");
    double geometry = par_file.get_parameter<double>("geometry");


    //double r_disc = par_file.get_parameter<double>("r_disc");

    if (par_args.key_exists("--source_h")) source[1] = par_args.get_parameter<double>("--source_h");

	const double r_isco = kerr_isco<double>(spin, +1);
	if(r_min < 0) r_min = r_isco;
    double dr = (logbin_r) ? exp(log(r_disc / r_min) / (Nr)) : (r_disc - r_min) / (Nr);

    long num_primary_rays = (((cosalphamax - cosalpha0)/dcosalpha) * ((betamax - beta0)/dbeta));

    long *disc_rays;
    double *disc_r, *disc_area, *disc_flux, *disc_emis, *disc_redshift, *disc_time;

    long disc_count = 0;

    disc_rays = new long[Nr];
    disc_r = new double[Nr];
    disc_area = new double[Nr];
    disc_flux = new double[Nr];
    disc_emis = new double[Nr];
    disc_redshift = new double[Nr];
    disc_time = new double[Nr];

    for(int ir=0; ir<Nr; ir++)
    {
        disc_r[ir] = (logbin_r) ? r_min * pow(dr, ir) : r_min + ir*dr;
        disc_area[ir] = integrate_disc_area(disc_r[ir], (logbin_r) ? disc_r[ir]*dr : dr, spin);

        disc_rays[ir] = 0;
        disc_flux[ir] = 0;
        disc_emis[ir] = 0;
        disc_redshift[ir] = 0;
        disc_time[ir] = 0;
    }

	PointSource<double> raytrace_source(source, V, spin, TOL, dcosalpha, dbeta, cosalpha0, cosalphamax, beta0, betamax);

    ZDestination<double>* my_destination = new ZDestination<double>(thetalim, r_disc);
    //DelayedFlaredDisc<double>* my_destination = new DelayedFlaredDisc<double>(thetalim, r_disc, r_angle_disc_dis);
    //AngledDiscsDestination<double> *my_destination = new AngledDiscsDestination<double>(thetalim, thetalim, r_angle_disc_dis);
    //TorusDiscDestination<double>* my_destination = new TorusDiscDestination<double>(r_torus, r_disc, r_isco);
    //InclPortionDiscDestination<double>* my_destination = new InclPortionDiscDestination<double>(thetalim, thetalim, r_angle_disc_dis);
    //EllipseDiscDestination<double>* my_destination = new EllipseDiscDestination<double>(r_disc, r_isco, major_axis, minor_axis);
    //SinDiscDestination<double>* my_destination = new SinDiscDestination<double>(r_disc);
    //ShakuraDiscDestination<double>* my_destination = new ShakuraDiscDestination<double>(efficiency, edd, r_isco);

    raytrace_source.redshift_start();
	raytrace_source.run_raytrace(my_destination, 1.1 * dist, r_disc);
	raytrace_source.range_phi();
	raytrace_source.redshift(my_destination, -1);

	for (int ray = 0; ray < raytrace_source.get_count(); ray++)
	{
		if (raytrace_source.rays[ray].steps > 0)
		{
		    double x, y, z;
		    cartesian(x, y, z, raytrace_source.rays[ray].r, raytrace_source.rays[ray].theta, raytrace_source.rays[ray].phi, spin);
            //if(z < 1E-2 && raytrace_source.rays[ray].redshift > 0 && raytrace_source.rays[ray].r >= r_isco)
            if(raytrace_source.rays[ray].redshift > 0 && raytrace_source.rays[ray].r >= r_isco)
            {
                int ir = (logbin_r) ? static_cast<int>( log(raytrace_source.rays[ray].r / r_min) / log(dr)) : static_cast<int>((raytrace_source.rays[ray].r - r_min) / dr);

                if(ir >= 0 && ir < Nr)
                {
                    ++disc_rays[ir];

                    // primary flux to work out returning radiation normalisation
                    // fraction of rays from the primary source hitting this part of the disc
                    // multiplied by the redshift to obtain the photon arrival rate in the rest frame
                    disc_flux[ir] += 1/(num_primary_rays * pow(raytrace_source.rays[ray].redshift,1));

                    // emissivity in the rest frame of the disc material
                    disc_emis[ir] += 1/pow(raytrace_source.rays[ray].redshift,gamma);

                    disc_redshift[ir] += raytrace_source.rays[ray].redshift;
                    disc_time[ir] += raytrace_source.rays[ray].t;
                }

                ++disc_count;
            }
		}
	}

	for(int ir=0; ir<Nr; ir++)
    {
        disc_redshift[ir] /= disc_rays[ir];
        disc_time[ir] /= disc_rays[ir];
        disc_flux[ir] /= disc_area[ir];
        disc_emis[ir] /= disc_area[ir];
    }

	TextOutput outfile((char*)out_filename.c_str());
	for(int ir=0; ir<Nr; ir++)
    {
	    outfile << disc_r[ir]
	            << disc_area[ir]
	            << disc_rays[ir]
	            << disc_flux[ir]
	            << disc_emis[ir]
	            << disc_redshift[ir]
	            << disc_time[ir]
	            << endl;
    }
	outfile.close();

	cout << "Done" << endl;

	return 0;
}
