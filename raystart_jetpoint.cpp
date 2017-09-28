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

#include "common/jetpointsource.h"

void read_par(const char* filename, char* outfile, double* source, double& V, double& a, double& dcosalpha, double& dbeta, double& cosalpha0, double& cosalphamax, double& beta0, double& betamax);

int main( )
{
	// parameter configuration file
	const char par_filename[] = "../par/raystart_jetpoint.par";
	char out_filename[128];

	double source[4];
	double V, spin;

	double cosalpha0;
	double cosalphamax;
	double dcosalpha;

	double beta0;
	double betamax;
	double dbeta;

	double *t, *r, *theta, *phi, *redshift;
	double *pt, *pr, *ptheta, *pphi;
	int *steps;


	JetPointSource<double> *RaytraceSource;


	// read in parameters from file
	read_par(par_filename, out_filename, source, V, spin, dcosalpha, dbeta, cosalpha0, cosalphamax, beta0, betamax);

	cout << "*****" << endl;
	cout << "Source r = [" << source[0] << " , " << source[1] << " , " << source[2] << " , " << source[3] << "] mu" << endl;
	cout << "Source radial velocity V = " << V << " mu*c" << endl;
	cout << "Spin a = " << spin << endl;
	cout << "*****" << endl << endl;

	RaytraceSource = new JetPointSource<double>( source, V, spin, TOL, dcosalpha, dbeta, cosalpha0, cosalphamax );
	RaytraceSource->CalculateMomentum( );
	RaytraceSource->MapResults(steps, t, r, theta, phi, redshift);
	RaytraceSource->MapMomentum(pt, pr, ptheta, pphi);

	TextOutput outfile(out_filename);
	for(int i=0; i<RaytraceSource->GetCount(); i++)
	{
		if(steps[i] < 0) continue;

		outfile << t[i] << r[i] << theta[i] << phi[i]
		        << pt[i] << pr[i] << ptheta[i] << pphi[i]
		        << endl;
	}


	outfile.Close();

	delete RaytraceSource;

	cout << "Done" << endl;

	return 0;
}

void read_par(const char* filename, char* outfile, double* source, double& V, double& a, double& dcosalpha, double& dbeta, double& cosalpha0, double& cosalphamax, double& beta0, double& betamax)
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
      break;
    }

  file_par.close();
}


