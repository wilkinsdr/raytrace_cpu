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
    double rmax = par_file.get_parameter<double>("rmax", 1000);
	int show_progress = par_file.get_parameter<int>("show_progress", 1);

	double r_isco = KerrISCO<double>(spin, +1);

	HealpixPointSource<double> raytrace_source(source, V, spin, order);

	raytrace_source.RedshiftStart();
	raytrace_source.RunRaytrace(rmax, M_PI_2, show_progress);
	raytrace_source.Redshift(-1);

    raytrace_source.MapResults(steps, t, r, theta, phi, redshift);

    TextOutput outfile((char*)out_filename.c_str());
    for(int pix=0; pix<raytrace_source.GetNumPix(); pix++)
    {
        double pix_r = r[5*pix+4];
        double pix_theta = theta[5*pix+4];
        if(pix_r > r_isco || pix_theta > (M_PI_2 - 1E-2))
            outfile << pix << r[5*pix+4] << redshift[5*pix+4] << endl;
        else
            outfile << pix << 0 << 0 << endl;
    }
    outfile.Close();

	return 0;
}
