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

int main( )
{
	// parameter configuration file
	const char par_filename[] = "../par/trace_rays.par";
	char out_filename[128];

	double source[4];
	double V, spin;

	double cosalpha0;
	double cosalphamax;
	double dcosalpha;

	double beta0;
	double betamax;
	double dbeta;

	double r_max, theta_max, r_isco;

	int write_step;

	PointSource<double> *RaytraceSource;


	// read in parameters from file
	read_par(par_filename, out_filename, source, V, spin, dcosalpha, dbeta, cosalpha0, cosalphamax, beta0, betamax, r_max, theta_max, write_step);

	cout << "*****" << endl;
	cout << "Source r = [" << source[0] << " , " << source[1] << " , " << source[2] << " , " << source[3] << "] mu" << endl;
	cout << "Source angular velocity V = " << V << " mu*c" << endl;
	cout << "Spin a = " << spin << endl;
	cout << "*****" << endl << endl;

	TextOutput outfile(out_filename);

	RaytraceSource = new PointSource<double>( source, V, spin, TOL, dcosalpha, dbeta, cosalpha0, cosalphamax );
	RaytraceSource->RunRaytrace( r_max, theta_max, &outfile, write_step );

	delete RaytraceSource;

	outfile.Close();

	cout << "Done" << endl;

	return 0;
}

void read_par(const char* filename, char* outfile, double* source, double& V, double& a, double& dcosalpha, double& dbeta, double& cosalpha0, double& cosalphamax, double& beta0, double& betamax, double& r_max, double& theta_max, int& write_step)
{
  //
  // read in program parameters from file
  //
  // arguments:
  //   filename    char*     name of file to read parameters from
  //   (others)              store values of parameters for main()
  //
  ifstream file_par;

  file_par.open(filename);

  while( !file_par.eof() )
    {
      file_par >> outfile;
      file_par >> source[0] >> source[1] >> source[2] >> source[3];
      file_par >> V;
      file_par >> a;
      file_par >> dcosalpha >> dbeta;
      file_par >> cosalpha0 >> cosalphamax;
      file_par >> beta0 >> betamax;
      file_par >> r_max >> theta_max;
      file_par >> write_step;
      break;
    }

  file_par.close();
}


