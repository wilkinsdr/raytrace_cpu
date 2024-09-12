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

	cout << "Setting up point source with " << Raytracer<T>::nRays << " rays" << endl;
    init_pointsource(pos, dcosalpha, dbeta, cosalpha0, cosalphamax, beta0, betamax);
}

//template <typename T>
//PointSource<T>::~PointSource()
//{
//
//}

template <typename T>
void PointSource<T>::init_pointsource(T* pos, T dcosalpha, T dbeta, T cosalpha0, T cosalphamax, T beta0, T betamax )
{
	for(int i=0; i<n_cosalpha; i++)
		for(int j=0; j<n_beta; j++)
		{
			int ix = i*n_beta + j;

			const T cosalpha = cosalpha0 + i*dcosalpha;
			const T beta = beta0 + j*dbeta;

			if(cosalpha >= cosalphamax || beta >= betamax)
			{
				Raytracer<T>::rays[ix].steps = -1;
				continue;
			}

			const T alpha = acos(cosalpha);

            Raytracer<T>::rays[ix].alpha = cosalpha;
            Raytracer<T>::rays[ix].beta = beta;

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
			Raytracer<T>::calculate_constants(ix, alpha, beta, velocity, energy);
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
void PointSource<T>::redshift(RayDestination<T>* destination, T V)
{
	//
	// Call the redshift_start function of the base class using the angular velocity for a circular orbit at the ray's end point
	// for rays incident on the accretion disc
	//
    Raytracer<T>::redshift(destination, V);
}

template class PointSource<double>;
template class PointSource<float>;
