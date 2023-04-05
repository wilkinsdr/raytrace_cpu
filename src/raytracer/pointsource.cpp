/*
 * pointsource.cpp
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#include "pointsource.h"

template <typename T>
PointSource<T>::PointSource( T* pos, T V, T spin, T tol, T dcosalpha, T dbeta, T cosalpha0, T cosalphamax, T beta0, T betamax, T E )
	: Raytracer<T>( (((cosalphamax - cosalpha0) / dcosalpha) + 1) * (((betamax - beta0) / dbeta) + 1) , spin , tol ),
	  velocity(V),
	  energy(E)
{
	n_cosalpha = ((cosalphamax - cosalpha0) / dcosalpha) + 1;
	n_beta = ((betamax - beta0) / dbeta) + 1;

	m_cosalpha = new T[Raytracer<T>::nRays];
	m_beta = new T[Raytracer<T>::nRays];

	cout << "Setting up point source with " << Raytracer<T>::nRays << " rays" << endl;
    init_pointsource(pos, dcosalpha, dbeta, cosalpha0, cosalphamax, beta0, betamax);
}

template <typename T>
PointSource<T>::~PointSource()
{
	delete[] m_cosalpha;
	delete[] m_beta;
}

template <typename T>
void PointSource<T>::init_pointsource(T* pos, T dcosalpha, T dbeta, T cosalpha0, T cosalphamax, T beta0, T betamax )
{
	for(int i=0; i<n_cosalpha; i++)
		for(int j=0; j<n_beta; j++)
		{
			int ix = i*n_beta + j;

			m_cosalpha[ix] = cosalpha0 + i*dcosalpha;
			m_beta[ix] = beta0 + j*dbeta;

			if(m_cosalpha[ix] >= cosalphamax || m_beta[ix] >= betamax)
			{
				Raytracer<T>::rays[ix].steps = -1;
				continue;
			}

			const T alpha = acos(m_cosalpha[ix]);

//			Raytracer<T>::m_rdot_sign[ix] = (alpha < M_PI/2) ? 1 : -1;
//			Raytracer<T>::m_thetadot_sign[ix] = (abs(m_beta[ix]) < M_PI/2) ? 1 : -1;

			// initialise position of photon
			Raytracer<T>::rays[ix].t = pos[0];
			Raytracer<T>::rays[ix].r = pos[1];
			Raytracer<T>::rays[ix].theta = pos[2];
			Raytracer<T>::rays[ix].phi = pos[3];;
			Raytracer<T>::rays[ix].pt = 0;
			Raytracer<T>::rays[ix].pr = 0;
			Raytracer<T>::rays[ix].ptheta = 0;
			Raytracer<T>::rays[ix].pphi = 0;
			Raytracer<T>::rays[ix].steps = 0;
			// calculate constants of motion
			Raytracer<T>::calculate_constants(ix, alpha, m_beta[ix], velocity, energy);
		}
}

template <typename T>
void PointSource<T>::redshift_start( )
{
	//
	// Call the redshift_start function of the base class using the source's angular velocity
	//
    Raytracer<T>::redshift_start(velocity);
}

template <typename T>
void PointSource<T>::redshift(T V )
{
	//
	// Call the redshift_start function of the base class using the angular velocity for a circular orbit at the ray's end point
	// for rays incident on the accretion disc
	//
    Raytracer<T>::redshift(V);
}

template class PointSource<double>;
template class PointSource<float>;
