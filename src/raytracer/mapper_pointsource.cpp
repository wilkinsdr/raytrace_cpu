/*
 * pointsource.cpp
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#include "mapper_pointsource.h"

template <typename T>
Mapper_PointSource<T>::Mapper_PointSource( T* pos, T V, T spin, T dcosalpha, T dbeta, T init_r0, T init_rmax, int init_Nr, int init_Ntheta, int init_Nphi, bool init_logbin_r, T init_thetamax, T cosalpha0, T cosalphamax, T beta0, T betamax, T tol, T E)
	: Mapper<T>( (((cosalphamax - cosalpha0) / dcosalpha) + 1) * (((betamax - beta0) / dbeta) + 1) , spin , init_r0, init_rmax, init_Nr, init_Ntheta, init_Nphi, init_logbin_r, init_thetamax, tol ),
	  velocity(V),
	  energy(E)
{
	n_cosalpha = ((cosalphamax - cosalpha0) / dcosalpha) + 1;
	n_beta = ((betamax - beta0) / dbeta) + 1;

	m_cosalpha = new T[Raytracer<T>::nRays];
	m_beta = new T[Raytracer<T>::nRays];

	cout << "Setting up point source with " << Raytracer<T>::nRays << " rays" << endl;
	InitPointSource( pos, dcosalpha, dbeta, cosalpha0, cosalphamax, beta0, betamax );
}

template <typename T>
Mapper_PointSource<T>::~Mapper_PointSource()
{
	delete[] m_cosalpha;
	delete[] m_beta;
}

template <typename T>
void Mapper_PointSource<T>::InitPointSource( T* pos, T dcosalpha, T dbeta, T cosalpha0, T cosalphamax, T beta0, T betamax )
{
	for(int i=0; i<n_cosalpha; i++)
		for(int j=0; j<n_beta; j++)
		{
			int ix = i*n_beta + j;

			m_cosalpha[ix] = cosalpha0 + i*dcosalpha;
			m_beta[ix] = beta0 + j*dbeta;

			if(m_cosalpha[ix] >= cosalphamax || m_beta[ix] >= betamax)
			{
				Raytracer<T>::m_steps[ix] = -1;
				return;
			}

			const T alpha = acos(m_cosalpha[ix]);

			Raytracer<T>::m_rdot_sign[ix] = (alpha < M_PI/2) ? 1 : -1;
			Raytracer<T>::m_thetadot_sign[ix] = (abs(m_beta[ix]) < M_PI/2) ? 1 : -1;

			// initialise position of photon
			Raytracer<T>::m_t[ix] = pos[0];
			Raytracer<T>::m_r[ix] = pos[1];
			Raytracer<T>::m_theta[ix] = pos[2];
			Raytracer<T>::m_phi[ix] = pos[3];;
			Raytracer<T>::m_pt[ix] = 0;
			Raytracer<T>::m_pr[ix] = 0;
			Raytracer<T>::m_ptheta[ix] = 0;
			Raytracer<T>::m_pphi[ix] = 0;
			Raytracer<T>::m_steps[ix] = 0;
			// calculate constants of motion
			Raytracer<T>::CalculateConstants(ix, alpha, m_beta[ix], velocity, energy);

		}
}

template <typename T>
void Mapper_PointSource<T>::RedshiftStart( )
{
	//
	// Call the RedshiftStart function of the base class using the source's angular velocity
	//
	Raytracer<T>::RedshiftStart( velocity );
}

template <typename T>
void Mapper_PointSource<T>::Redshift( T V )
{
	//
	// Call the RedshiftStart function of the base class using the angular velocity for a circular orbit at the ray's end point
	// for rays incident on the accretion disc
	//
	Raytracer<T>::Redshift( V );
}

template class Mapper_PointSource<double>;
template class Mapper_PointSource<float>;