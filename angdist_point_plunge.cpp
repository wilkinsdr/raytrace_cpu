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

#include "common/pointsource_vel.h"
#include "common/par_file.h"

int main(int argc, char** argv)
{
	// parameter configuration file
	char default_par_filename[] = "../par/angdist_point_plunge.par";
	char* par_filename;
	if(argc == 2)
		par_filename = argv[1];
	else
		par_filename = default_par_filename;

	double source[4];
	double V[4];

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
	int Ntheta;

	double *btheta, *count, *count_perarea;


	PointSourceVel<double> *RaytraceSource;

	ParameterFile par_file(par_filename);
	string out_filename = par_file.get_parameter<string>("outfile");
	double source_r = par_file.get_parameter<double>("source_r");
    double source_phi = par_file.get_parameter<double>("source_phi", 1.5707);
	double spin = par_file.get_parameter<double>("spin");
	double cosalpha0 = par_file.get_parameter<double>("cosalpha0", -0.9999);
	double cosalphamax = par_file.get_parameter<double>("cosalphamax", 0.9999);
	double dcosalpha = par_file.get_parameter<double>("dcosalpha");
	double beta0 = par_file.get_parameter<double>("beta0", -1*M_PI);
	double betamax = par_file.get_parameter<double>("betamax", M_PI);
	double dbeta = par_file.get_parameter<double>("dbeta");
	double danglebin = par_file.get_parameter<double>("danglebin");
	int axis = par_file.get_parameter<int>("axis");
	int axis2 = par_file.get_parameter<int>("axis2");

	if(axis < 1 || axis > 3 || abs(axis2) < 1 || abs(axis2) > 3)
	{
		cerr << "Parameter ERROR: Axis must be 1, 2 or 3" << endl;
		return 1;
	}

	double r_isco = KerrISCO(spin, +1);

    //const T rhosq = r*r + (spin*__cosf(disc_theta))*(spin*__cosf(disc_theta));
    const double delta = source_r*source_r - 2*source_r + spin*spin;
    //const T sigmasq = (r*r + spin*spin)*(r*r + spin*spin) - spin*spin*delta*__sinf(disc_theta)*__sinf(disc_theta);

    // constants of motion for the ISCO
    const double u = 1/r_isco;
    const double k = (1 - 2*u + spin*u*sqrtf(u)) / sqrtf(1 - 3*u + 2*spin*u*sqrtf(u));
    const double h = (1 + spin*spin*u*u - 2*spin*u*sqrtf(u)) / sqrtf( u*(1 - 3*u + 2*spin*u*sqrtf(u)) );

    // 4-velocity of a plunging orbit from the ISCO
    V[0] = (1/delta) * ( (source_r*source_r + spin*spin + 2*spin*spin/source_r)*k - 2*spin*h/source_r );
    V[1] = -1*sqrtf( k*k - 1 + 2/source_r + (spin*spin*(k*k - 1) - h*h)/(source_r*source_r) + 2*(h-spin*k)*(h-spin*k)/(source_r*source_r*source_r) );
    V[2] = 0;
    V[3] = (1/delta) * (2*spin*k/source_r + (1 - 2/source_r)*h);

    source[0] = 0;
    source[1] = source_r;
    source[2] = M_PI/2 - 0.001;
    source[3] = source_phi;
    

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

	RaytraceSource = new PointSourceVel<double>( source, V, spin, TOL, dcosalpha, dbeta, cosalpha0, cosalphamax );
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

