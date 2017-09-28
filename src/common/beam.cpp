/*
 * Beam.cpp
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#include "Beam.h"

template <typename T>
Beam<T>::Beam( T* pos, T V, T spin, T tol, T cosalpha0, T beta0, T dcosalpha, T dbeta,  T E )
	: Raytracer<T>( 4 , spin , tol ),
	  velocity(V),
	  energy(E)
{
	n_cosalpha = ((cosalphamax - cosalpha0) / dcosalpha) + 1;
	n_beta = ((betamax - beta0) / dbeta) + 1;

	m_cosalpha = new T[Raytracer<T>::nRays];
	m_beta = new T[Raytracer<T>::nRays];

	cout << "Setting up beam with " << Raytracer<T>::nRays << " rays" << endl;
	InitBeam( pos, dcosalpha, dbeta, cosalpha0, cosalphamax, beta0, betamax );
}

template <typename T>
Beam<T>::~Beam()
{
	delete[] m_cosalpha;
	delete[] m_beta;
}

template <typename T>
void Beam<T>::InitBeam( T* pos, T dcosalpha, T dbeta, T cosalpha0, T cosalphamax, T beta0, T betamax )
{
	for(int i=0; i<n_cosalpha; i++)
		for(int j=0; j<n_beta; j++)
		{
			int ix = i*n_beta + j;

			switch(i)
			{
				case 0:
					m_cosalpha[ix] = cosalpha0 + 0.5*dcosalpha;
					m_beta[ix] = beta0 + 0.5*dbeta;
					break;
				case 1:
					m_cosalpha[ix] = cosalpha0 - 0.5*dcosalpha;
					m_beta[ix] = beta0 + 0.5*dbeta;
					break;
				case 2:
					m_cosalpha[ix] = cosalpha0 - 0.5*dcosalpha;
					m_beta[ix] = beta0 - 0.5*dbeta;
					break;
				case 3:
					m_cosalpha[ix] = cosalpha0 + 0.5*dcosalpha;
					m_beta[ix] = beta0 - 0.5*dbeta;
					break;

			}
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
void Beam<T>::RedshiftStart( )
{
	//
	// Call the RedshiftStart function of the base class using the source's angular velocity
	//
	Raytracer<T>::RedshiftStart( velocity );
}

template <typename T>
void Beam<T>::Redshift( T V )
{
	//
	// Call the RedshiftStart function of the base class using the angular velocity for a circular orbit at the ray's end point
	// for rays incident on the accretion disc
	//
	Raytracer<T>::Redshift( V );
}

template class Beam<double>;
template class Beam<float>;
