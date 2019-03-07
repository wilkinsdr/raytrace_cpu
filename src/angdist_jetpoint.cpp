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

#include "raytracer/jetpointsource.h"

void read_par(const char* filename, char* outfile, double* source, double& V, double& a, double& dcosalpha, double& dbeta, double& cosalpha0, double& cosalphamax, double& beta0, double& betamax, double& dthetabin, int& axis, int& axis2);

int main( )
{
	// parameter configuration file
	const char par_filename[] = "../par/angdist_jetpoint.par";
	char out_filename[128];

	double source[4];
	double V, spin;

	double cosalpha0;
	double cosalphamax;
	double dcosalpha;

	double beta0;
	double betamax;
	double dbeta;

	double et[4];
	double e1[4];
	double e2[4];
	double e3[4];
	double g[4][4];
	double mink[4][4];

	double *t, *r, *theta, *phi, *redshift;
	double *pt, *pr, *ptheta, *pphi;
	int *steps;

	double anglebin0 = 0.0;
	double anglebinmax = 2*M_PI;
	double danglebin;
	int Ntheta;

	int axis, axis2;

	double *btheta, *count, *count_perarea;


	JetPointSource<double> *RaytraceSource;


	// read in parameters from file
	read_par(par_filename, out_filename, source, V, spin, dcosalpha, dbeta, cosalpha0, cosalphamax, beta0, betamax, danglebin, axis, axis2);

	if(axis < 1 || axis > 3 || abs(axis2) < 1 || abs(axis2) > 3)
	{
		cerr << "Parameter ERROR: Axis must be 1, 2 or 3" << endl;
		return 1;
	}

	cout << "*****" << endl;
	cout << "Source r = [" << source[0] << " , " << source[1] << " , " << source[2] << " , " << source[3] << "] mu" << endl;
	cout << "Source radial velocity V = " << V << " mu*c" << endl;
	cout << "Spin a = " << spin << endl;
	cout << "*****" << endl << endl;

	Ntheta = static_cast<int>((anglebinmax - anglebin0) / danglebin) + 1;
	btheta = new double[Ntheta];
	count = new double[Ntheta];
	count_perarea = new double[Ntheta];
	for(int i=0; i<Ntheta; i++)
	{
		btheta[i] = anglebin0 + i*danglebin;
		count[i] = 0;
	}

	// get the metric and the basis vectors at the location of the source for a stationary observer
	Metric(g, source, spin);
	Tetrad<double>(et, e1, e2, e3, source, 0., spin);
	// and the Minkowski metric
	Minkowski<double>(mink);

	RaytraceSource = new JetPointSource<double>( source, V, spin, TOL, dcosalpha, dbeta, cosalpha0, cosalphamax );
	RaytraceSource->CalculateMomentum( );
	RaytraceSource->MapResults(steps, t, r, theta, phi, redshift);
	RaytraceSource->MapMomentum(pt, pr, ptheta, pphi);

	for(int i=0; i<RaytraceSource->GetCount(); i++)
	{
		if(steps[i] < 0) continue;

		double p[] = { pt[i], pr[i], ptheta[i], pphi[i] };

		// transform the vector into the frame of a stationary observer at the source
		double pprime[] = { DotProduct(g, p, et), DotProduct(g, p, e1), DotProduct(g, p, e2), DotProduct(g, p, e3) };

//		cout << setw(15) << DotProduct(g, p, p) << setw(15) << DotProduct(mink, pprime, pprime) << endl;

		double cosangle = pprime[axis] / sqrt(pprime[1]*pprime[1] + pprime[2]*pprime[2] + pprime[3]*pprime[3]);
		double angle = acos(cosangle);

		if(axis2 > 0 && pprime[axis2] < 0) angle = 2*M_PI - angle;
		if(axis2 < 0 && pprime[-1*axis2] > 0) angle = 2*M_PI - angle;

		int anglebin = (angle - anglebin0) / danglebin;
		if(anglebin < 0 || anglebin >= Ntheta) continue;
		++count[anglebin];
	}

	double maxcount = 0;
	double maxcountperarea = 0;
	for(int i=0; i<Ntheta; i++)
	{
		count_perarea[i] = count[i]/(2*M_PI * danglebin * abs(sin(btheta[i] + 0.5*danglebin)));
		if(count[i] > maxcount) maxcount = count[i];
		if(count_perarea[i] > maxcountperarea && !isinf((double)count_perarea[i])) maxcountperarea = count_perarea[i];
	}

	TextOutput outfile(out_filename);
	for(int i=0; i<Ntheta; i++)
	{
		if( btheta[i] < 0.1 || btheta[i] > 2*M_PI-0.1 || (btheta[i] > M_PI-0.1 && btheta[i] < M_PI+0.1) ) continue;
		outfile << btheta[i] << count[i]/maxcount << count_perarea[i]/maxcountperarea << endl;
	}
	outfile.Close();

	delete[] btheta;
	delete[] count;

	delete RaytraceSource;

	cout << "Done" << endl;

	return 0;
}

void read_par(const char* filename, char* outfile, double* source, double& V, double& a, double& dcosalpha, double& dbeta, double& cosalpha0, double& cosalphamax, double& beta0, double& betamax, double& dthetabin, int& axis, int& axis2)
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
      file_par >> dthetabin;
      file_par >> axis >> axis2;
      break;
    }

  file_par.close();
}


