/*
 * pointsource.cpp
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#include "mapper_pointsource.h"


Mapper_PointSource::Mapper_PointSource( double* pos, double V, double spin, double dcosalpha, double dbeta, double init_r0, double init_rmax, int init_Nr, int init_Ntheta, int init_Nphi, bool init_logbin_r, double init_thetamax, double cosalpha0, double cosalphamax, double beta0, double betamax, double tol, double E)
	: Mapper( (((cosalphamax - cosalpha0) / dcosalpha) + 1) * (((betamax - beta0) / dbeta) + 1) , spin , init_r0, init_rmax, init_Nr, init_Ntheta, init_Nphi, init_logbin_r, init_thetamax, tol ),
	  velocity(V),
	  energy(E)
{
	n_cosalpha = ((cosalphamax - cosalpha0) / dcosalpha) + 1;
	n_beta = ((betamax - beta0) / dbeta) + 1;

	m_cosalpha = new double[Raytracer<double>::nRays];
	m_beta = new double[Raytracer<double>::nRays];

	cout << "Setting up point source with " << Raytracer<double>::nRays << " rays" << endl;
	InitPointSource( pos, dcosalpha, dbeta, cosalpha0, cosalphamax, beta0, betamax );
}


Mapper_PointSource::~Mapper_PointSource()
{
	delete[] m_cosalpha;
	delete[] m_beta;
}


void Mapper_PointSource::InitPointSource( double* pos, double dcosalpha, double dbeta, double cosalpha0, double cosalphamax, double beta0, double betamax )
{
	for(int i=0; i<n_cosalpha; i++)
		for(int j=0; j<n_beta; j++)
		{
			int ix = i*n_beta + j;

			m_cosalpha[ix] = cosalpha0 + i*dcosalpha;
			m_beta[ix] = beta0 + j*dbeta;

			if(m_cosalpha[ix] >= cosalphamax || m_beta[ix] >= betamax)
			{
				Raytracer<double>::m_steps[ix] = -1;
				return;
			}

			const double alpha = acos(m_cosalpha[ix]);

			Raytracer<double>::m_rdot_sign[ix] = (alpha < M_PI/2) ? 1 : -1;
			Raytracer<double>::m_thetadot_sign[ix] = (abs(m_beta[ix]) < M_PI/2) ? 1 : -1;

			// initialise position of photon
			Raytracer<double>::m_t[ix] = pos[0];
			Raytracer<double>::m_r[ix] = pos[1];
			Raytracer<double>::m_theta[ix] = pos[2];
			Raytracer<double>::m_phi[ix] = pos[3];;
			Raytracer<double>::m_pt[ix] = 0;
			Raytracer<double>::m_pr[ix] = 0;
			Raytracer<double>::m_ptheta[ix] = 0;
			Raytracer<double>::m_pphi[ix] = 0;
			Raytracer<double>::m_steps[ix] = 0;
			// calculate constants of motion
			Raytracer<double>::CalculateConstants(ix, alpha, m_beta[ix], velocity, energy);

		}
}


void Mapper_PointSource::RedshiftStart( )
{
	//
	// Call the RedshiftStart function of the base class using the source's angular velocity
	//
	Raytracer<double>::RedshiftStart( velocity );
}


void Mapper_PointSource::Redshift( double V )
{
	//
	// Call the RedshiftStart function of the base class using the angular velocity for a circular orbit at the ray's end point
	// for rays incident on the accretion disc
	//
	Raytracer<double>::Redshift( V );
}

