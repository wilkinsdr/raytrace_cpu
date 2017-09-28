/*
 * trace_rays.cpp
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

#include "common/pointsource.h"
#include "common/par_file.h"

void read_par(const char* filename, char* outfile, double* source, double& V, double& a, double& dcosalpha, double& dbeta, double& cosalpha0, double& cosalphamax, double& beta0, double& betamax, double& r_max, double& theta_max, int& write_step);

int main(int argc, char** argv)
{
	// parameter configuration file
	char default_par_filename[] = "../par/pointsource_sky_discfrac.par";
	char* par_filename;
	if(argc == 2)
		par_filename = argv[1];
	else
		par_filename = default_par_filename;

	double source[4];
	double h0, hmax, dh;
	int Nh;

	double V, spin;

	double cosalpha0;
	double cosalphamax;
	double dcosalpha;

	double beta0;
	double betamax;
	double dbeta;

	double r_disc, r_max, theta_disc, r_isco;

	double *t, *r, *theta, *phi, *redshift;
	int* steps;
	double domega, omega_disc, omega_total;

	ParameterFile par_file(par_filename);
	out_filename = par_file.get_parameter<string>("outfile");
	h0 = par_file.get_parameter<double>("h0");
	hmax = par_file.get_parameter<double>("hmax");
	Nh = par_file.get_parameter<int>("Nh");
	logbin_h = par_file.get_parameter<bool>("logbin_h");
	V = par_file.get_parameter<double>("V");
	spin = par_file.get_parameter<double>("spin");
	cosalpha0 = par_file.get_parameter<double>("cosalpha0", -0.9999);
	cosalphamax = par_file.get_parameter<double>("cosalphamax", 0.9999);
	dcosalpha = par_file.get_parameter<double>("dcosalpha");
	beta0 = par_file.get_parameter<double>("beta0", -1*M_PI);
	betamax = par_file.get_parameter<double>("betamax", M_PI);
	dbeta = par_file.get_parameter<double>("dbeta");
	r_disc = par_file.get_parameter<double>("r_disc", 400);
	theta_disc = par_file.get_parameter<double>("theta_disc", M_PI/2);
	r_max = par_file.get_parameter<double>("r_max", 1000);

	PointSource<double> *RaytraceSource;

	cout << "*****" << endl;
	cout << "Source r = [" << source[0] << " , " << source[1] << " , " << source[2] << " , " << source[3] << "] mu" << endl;
	cout << "Source angular velocity V = " << V << " mu*c" << endl;
	cout << "Spin a = " << spin << endl;
	cout << "*****" << endl << endl;

	r_isco = KerrISCO(spin, +1);

	source[0] = 0;
	source[1] = 0;
	source[2] = 1E-3;
	source[3] = 4.712;

	TextOutput outfile((const char*)out_filename.c_str());

	for(int ih = 0; ih < Nh; ih++)
	{
		source[1] = (logbin_h) ? h0 * pow(dh, ih) : r0 + ih * dh;

		cout << "Source at h = " << source[1] << endl;


		RaytraceSource = new PointSource<double>(source, V, spin, TOL, dcosalpha, dbeta, cosalpha0, cosalphamax);
		RaytraceSource->RunRaytrace(r_max, theta_disc);
		RaytraceSource->MapResults(steps, t, r, theta, phi, redshift);

		omega_disc = 0;
		omega_total = 0;

		for (int i = 0; i < RaytraceSource->GetCount(); i++) {
			domega = dcosalpha * dbeta;

			omega_total += domega;

			if (steps[i] > 0) {
				if (theta[i] >= theta_disc && r[i] >= r_isco) {
					omega_disc += domega;
				}
			}
		}

		delete RaytraceSource;

		outfile << source[1] << omega_disc << omega_disc / omega_total << endl;
	}

	outfile.Close();

	cout << "Done" << endl;

	return 0;
}

