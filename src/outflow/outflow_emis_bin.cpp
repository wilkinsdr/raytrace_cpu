/*
 * trace_rays.cpp
 *
 *  Created on: 7 Mar 2019
 *      Author: drw
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

#include "mapper/mapper_imageplane.h"
#include "include/par_file.h"

bool pointsource_stop(double t, double r, double theta, double phi, double* args)
{
	const double a = args[0];
	const double h = args[1];
	const double rad = args[2];

	const double x = sqrt(r*r + a*a)*sin(theta)*cos(phi);
	const double y = sqrt(r*r + a*a)*sin(theta)*sin(phi);
	const double z = r*cos(theta);

	if( (x*x + y*y + (z-h)*(z-h)) < rad ) return true;

	//if(r < 2) return true;

	return false;
}

int main(int argc, char** argv)
{
	// parameter configuration file
	char default_par_filename[] = "../par/outflow_emis_bin.par";
	char* par_filename;
	if(argc == 2)
		par_filename = argv[1];
	else
		par_filename = default_par_filename;

	Mapper_ImagePlane *emis_mapper;

	ParameterFile par_file(par_filename);
	string out_filename = par_file.get_parameter<string>("outfile");
	double dist = par_file.get_parameter<double>("dist");
	double incl = par_file.get_parameter<double>("incl");
	double plane_phi0 = par_file.get_parameter<double>("plane_phi0", 0);
	double spin = par_file.get_parameter<double>("spin");
	double x0 = par_file.get_parameter<double>("x0");
	double xmax = par_file.get_parameter<double>("xmax");
	int Nx = par_file.get_parameter<int>("Nx");
	double y0 = par_file.get_parameter<double>("y0");
	double ymax = par_file.get_parameter<double>("ymax");
	int Ny = par_file.get_parameter<int>("Ny");
	double en0 = par_file.get_parameter<double>("en0");
	double enmax = par_file.get_parameter<double>("enmax");
	int Nen = par_file.get_parameter<int>("Nen");
	bool logbin_en = par_file.get_parameter<bool>("logbin_en", false);
	double outflow_rin = par_file.get_parameter<double>("outflow_rin");
	double outflow_rout = par_file.get_parameter<double>("outflow_rout");
    double outflow_thetamin = par_file.get_parameter<double>("outflow_thetamin", 0);
    double outflow_thetamax = par_file.get_parameter<double>("outflow_thetamax", M_PI_2);
    int Nr = par_file.get_parameter<int>("Nr", 50);
    bool logbin_r = par_file.get_parameter<bool>("logbin_r", false);
    int Ntheta = par_file.get_parameter<int>("Ntheta", 50);
    int Nphi = par_file.get_parameter<int>("Nphi", 50);
	double source_vel = par_file.get_parameter<double>("source_vel");
	int source_motion = par_file.get_parameter<int>("source_motion", 1);
	double tol = par_file.get_parameter<double>("tol", TOL);
    int show_progress = par_file.get_parameter<double>("show_progress", 1);

	double dx = (xmax - x0) / (Nx - 1);
	double dy = (ymax - y0) / (Ny - 1);

	cout << "*****" << endl;
	cout << "Image plane at  d = " << dist << " , incl = " << incl << endl;
	cout << "Spin a = " << spin << endl;
	cout << "*****" << endl << endl;

	emis_mapper = new Mapper_ImagePlane(dist, incl, x0, xmax, dx, y0, ymax, dy, spin, outflow_rin, outflow_rout, Nr, Ntheta, Nphi, logbin_r, M_PI_2, tol, plane_phi0);
	emis_mapper->set_motion(source_vel, source_motion);
	emis_mapper->redshift_start();
	emis_mapper->run_map( 1.5*dist, show_progress );

	const double den = (logbin_en) ? exp(log(enmax/en0)/Nen) : (enmax - en0)/Nen;

	double *line;
	line = new double[Nen];
    for(int ien=0; ien<Nen; ien++)
        line[ien] = 0;

	const int ir_start = emis_mapper->r_index(outflow_rin);
    const int ir_end = emis_mapper->r_index(outflow_rout);
    const int itheta_start = emis_mapper->theta_index(outflow_thetamin);
    const int itheta_end = emis_mapper->theta_index(outflow_thetamax);

    cout << "Computing wind emission..." << endl;
    for(int ir=ir_start; ir<=ir_end; ir++)
    {
        const double r = emis_mapper->bin_r(ir);
        const double bin_dr = (emis_mapper->logbin_r) ? emis_mapper->bin_dr * r : emis_mapper->bin_dr;
        const double flux = 1. / (4. * M_PI * r * r);
        for(int itheta=itheta_start; itheta <= itheta_end; itheta++)
        {
            const double theta = emis_mapper->bin_theta(itheta);
            const double area = sin(theta) * emis_mapper->bin_dtheta * emis_mapper->bin_dphi;
            cout << "area = " << area << endl;
            const double vol = area * bin_dr;
            cout << "vol = " << vol << endl;
            const double emissivity = flux * area * vol;
            cout << "emissivity = " << emissivity << endl;
            for(int iphi=0; iphi<emis_mapper->Nphi; iphi++)
            {
                if(emis_mapper->map_Nrays->elem(ir,itheta,iphi) > 0)
                {
                    const double enshift = 1. / emis_mapper->map_redshift->elem(ir,itheta,iphi);
                    const int ien = (logbin_en) ? static_cast<int>( log(enshift / en0) / log(den)) : static_cast<int>(
                            (enshift - en0) / den);
                    if (ien >= 0 && ien < Nen)
                        line[ien] += emissivity * pow(enshift, 3);
                }
            }
        }
    }

	TextOutput outfile((const char*)out_filename.c_str());
	for(int ien=0; ien<Nen; ien++)
	{
	    const double en = (logbin_en) ? en0 * pow(den, ien) : en0 + den * ien;
		outfile << en << line[ien] << endl;
	}
    outfile.close();

	cout << "Done" << endl;

	return 0;
}
