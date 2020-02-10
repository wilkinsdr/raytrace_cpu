/*
 * pointsource.cpp
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#include "jetpointsource.h"

template <typename T>
JetPointSource<T>::JetPointSource( T* pos, T V, T spin, T tol, T dcosalpha, T dbeta, T cosalpha0, T cosalphamax, T beta0, T betamax, T E )
	: Raytracer<T>( (((cosalphamax - cosalpha0) / dcosalpha) + 1) * (((betamax - beta0) / dbeta) + 1) , spin , tol ),
	  velocity(V),
	  energy(E)
{
	n_cosalpha = ((cosalphamax - cosalpha0) / dcosalpha) + 1;
	n_beta = ((betamax - beta0) / dbeta) + 1;

	m_cosalpha = new T[Raytracer<T>::nRays];
	m_beta = new T[Raytracer<T>::nRays];

	cout << "Setting up point source with " << Raytracer<T>::nRays << " rays" << endl;
	InitJetPointSource( pos, dcosalpha, dbeta, cosalpha0, cosalphamax, beta0, betamax );
}

template <typename T>
JetPointSource<T>::~JetPointSource()
{
	delete[] m_cosalpha;
	delete[] m_beta;
}

template <typename T>
void JetPointSource<T>::init_jet_pointsource(T* pos, T dcosalpha, T dbeta, T cosalpha0, T cosalphamax, T beta0, T betamax )
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

			// rdot_sign and theta_dot sign are sett in CalculateConstants_RadialMotion

			//cerr << velocity << ' ' << m_cosalpha[ix] << ' ' << Raytracer<T>::m_rdot_sign[ix] << endl;

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
			CalculateConstants_RadialMotion(ix, alpha, m_beta[ix], velocity, energy);

		}
}

template <typename T>
void JetPointSource<T>::redshift_start(  )
{
	//
	// Calculates the initial energy of each ray at emission for use in redshift calculations.
	// Emitter is orbitting the black hole rotation axis at angular velocity V, or if V = -1, the
	// angular velocity is calculated for a circular orbit at the current r co-ordinate of the ray.
	//
	// If ray is being propagated backwards in time (for image planes), the spin parameter is revered
	// (back to its true value as it will already have been reversed for the propagation) and the
	// spatial components of the ray's 4-momentum will be reversed so the ray is travelling in the
	// correct direction wrt the emitting material.
	//
	// This function should be called prior to running the raytrace if redshifts are required.
	//
	// Calls the GPURedshiftStart kernel to perform calciulation on the GPU.
	//
	// Arguments:
	//	V		T		Angular velocity of emitter (set to -1 for circular orbit at current r of ray)
	//	reverse	bool	Whether the ray is being propagated backwards in time (for image planes, default value = false)
	//
	cout << "Calculating initial energies" << endl;

	for(int ray = 0; ray < Raytracer<T>::nRays; ray++)
	{
		T p[4];
		const T V = velocity;
		const T a = Raytracer<T>::spin;

		const T t = Raytracer<T>::m_t[ray];
		const T r = Raytracer<T>::m_r[ray];
		const T theta = Raytracer<T>::m_theta[ray];
		const T phi = Raytracer<T>::m_phi[ray];
		const T k = Raytracer<T>::m_k[ray];
		const T h = Raytracer<T>::m_h[ray];
		const T Q = Raytracer<T>::m_Q[ray];
		const int rdot_sign = Raytracer<T>::m_rdot_sign[ray];
		const int thetadot_sign = Raytracer<T>::m_thetadot_sign[ray];

		// metric coefficients
		const T rhosq = r*r + (a*cos(theta))*(a*cos(theta));
		const T delta = r*r - 2*r + a*a;
		const T sigmasq = (r*r + a*a)*(r*r + a*a) - a*a*delta*sin(theta)*sin(theta);

		const T e2nu = rhosq * delta / sigmasq;
		const T e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
		const T omega = 2*a*r / sigmasq;

		T g[16];
		for(int i=0; i<16; i++)
			g[i] = 0;

		// g[i][j] -> g[i*4 + j]
		g[0*4 + 0] = e2nu - omega*omega*e2psi;
		g[0*4 + 3] = omega*e2psi;
		g[3*4 + 0] = g[0*4 + 3];
		g[1*4 + 1] = -rhosq/delta;
		g[2*4 + 2] = -rhosq;
		g[3*4 + 3] = -e2psi;

		// timelike basis vector
		const T et[] = { 1 / sqrtf(g[0] + g[5]*V*V), V / sqrtf(g[0] + g[5]*V*V), 0, 0 };

		// photon momentum
		MomentumFromConsts<T>(p[0], p[1], p[2], p[3], k, h, Q, rdot_sign, thetadot_sign, r, theta, phi, a);

		// evaluate dot product to get energy
		Raytracer<T>::m_emit[ray] = 0;
		for(int i=0; i<4; i++)
			for(int j=0; j<4; j++)
				Raytracer<T>::m_emit[ray] += g[i*4 + j] * et[i]* p[j];
	}
}

template <typename T>
void JetPointSource<T>::redshift(T V )
{
	//
	// Call the RedshiftStart function of the base class using the angular velocity for a circular orbit at the ray's end point
	// for rays incident on the accretion disc
	//
	Raytracer<T>::Redshift( V );
}

template <typename T>
inline void JetPointSource<T>::CalculateConstants_RadialMotion(int ray, T alpha, T beta, T V, T E)
{
	//
	// Compute the constants of motion for a ray emitted at polar angles alpha and beta in the frame
	// of a source at (t,r,theta,phi) moving radially at dr/dt = V
	//
	const T a = Raytracer<T>::spin;
	const T t = Raytracer<T>::m_t[ray];
	const T r = Raytracer<T>::m_r[ray];
	const T theta = Raytracer<T>::m_theta[ray];
	const T phi = Raytracer<T>::m_phi[ray];

	T k, h, Q;

	const T rhosq = r*r + (a*cos(theta))*(a*cos(theta));
	const T delta = r*r - 2*r + a*a;
	const T sigmasq = (r*r + a*a)*(r*r + a*a) - a*a*delta*sin(theta)*sin(theta);

	// metric coefficients
	const T e2nu = rhosq * delta / sigmasq;
	const T e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
	const T omega = 2*a*r / sigmasq;
	const T g00 = e2nu - omega*omega*e2psi;
	const T g03 = omega*e2psi;
	const T g11 = -rhosq/delta;
	const T g22 = -rhosq;
	const T g33 = -e2psi;

	// tetrad basis vector components (they're right-handed!)
	const T et0 = 1 / sqrt(g00 + g11*V*V);
	const T et1 = V * et0;
	//
	const T e12 = 1/sqrt(rhosq);
	//
	const T e23 = sqrt( g00 / (g03*g03 - g00*g33) );
	const T e20 = -(g03/g00) * e23;
	//
	const T e30 = sqrt( -g11 / g00 ) * V / sqrt( g00 + g11*V*V );
	const T e31 = sqrt( -g00 / g11 ) / sqrt( g00 + g11*V*V );

	// photon 4-momentum in source frame
	const T rdotprime[] = { E, E*sin(alpha)*cos(beta), E*sin(alpha)*sin(beta), E*cos(alpha) };

	const T tdot = rdotprime[0]*et0 + rdotprime[1]*e20 + rdotprime[3]*e30;
	const T phidot = rdotprime[1]*e23;
	const T rdot = rdotprime[0]*et1 + rdotprime[3]*e31;
	const T thetadot = rdotprime[2]*e12;

	// find the corresponding values of k, h and Q using the geodesic equations
	k = (1 - 2*r/rhosq)*tdot + (2*a*r*sin(theta)*sin(theta)/rhosq)*phidot;

	h = phidot * ( (r*r + a*a)*(r*r + a*a*cos(theta)*cos(theta) - 2*r)*sin(theta)*sin(theta) + 2*a*a*r*sin(theta)*sin(theta)*sin(theta)*sin(theta) );
	h = h - 2*a*r*k*sin(theta)*sin(theta);
	h = h / ( r*r + a*a*cos(theta)*cos(theta) - 2*r );

	Q = rhosq*rhosq*thetadot*thetadot - (a*k*cos(theta) + h/tan(theta))*(a*k*cos(theta) - h/tan(theta));

	Raytracer<T>::m_k[ray] = k;
	Raytracer<T>::m_h[ray] = h;
	Raytracer<T>::m_Q[ray] = Q;

	Raytracer<T>::m_rdot_sign[ray] = (rdot > 0) ? 1 : -1;
	Raytracer<T>::m_thetadot_sign[ray] = (thetadot > 0) ? 1 : -1;

//	// apparent velocity of source to a stationary observer at the same location
//	double vel = velocity * sqrt( -g11 / g00 );
//	// use the special relativistic abberation between source and stationary observer at the same location
//	// to work out the sign for rdot
//	// (work out what fraction dr/dt is of dr/dt for photons)
//	Raytracer<T>::m_rdot_sign[ray] = (cos(alpha) > -1*vel) ? 1 : -1;
//	Raytracer<T>::m_thetadot_sign[ray] = (abs(beta) < M_PI/2) ? 1 : -1;

}

template class JetPointSource<double>;
template class JetPointSource<float>;
