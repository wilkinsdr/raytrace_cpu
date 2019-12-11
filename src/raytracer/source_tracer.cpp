/*
 * raytracer.cu
 *
 *  Created on: 4 Sep 2013
 *      Author: drw
 */

#include "raytracer.h"

template <typename T>
Raytracer<T>::Raytracer( int num_rays, float spin_par, float toler, T en0, T enmax, int Nen, bool logbin_en )
	: Raytracer(num_rays, spin_par, toler)
{
	emis = new T*[num_rays];
	absorb = new T*[num_rays];

	for(int i=0; i<num_rays; i++)
	{
		emis[i] = new T[Nen];
		absorb[i] = new T[Nen];
		for(int j=0; j<Nen; j++)
		{
			emis[i][j] = 0;
			absorb[i][j] = 0;
		}
	}
}

template <typename T>
Raytracer<T>::~Raytracer( )
{

}

template <typename T>
void Raytracer<T>::RunRaytrace( T r_max, T theta_max, TextOutput* outfile, int write_step, T write_rmax, T write_rmin, bool write_cartesian)
{
	//
	// Runs the ray tracing algorithm once the rays have been set up.
	// Integration of geodesic equation for each ray proceeds until it reaches limiting radius or is lost through the event horizon.
	// Alternatively, Integration stop when the maximum number of steps (defined by STEPLIM in the header file) is reached to prevent
	// incredibly long loops of very small steps.
	//
	// To work around thread time limits on devices also running an X server, each kernel execution
	// only integrates for the number of steps defined in THREAD_STEPLIM. An unfinished flag is then
	// set and the kernel is executed repeatedly until each ray has finished.
	//
	// Calls the GPURaytrace kernel to perform calculation on the GPU.
	//
	// Arguments:
	//	r_max		T		Limiting outer radius for propagation in Rg (default value = 1000)
	//
	cout << "Running raytracer..." << endl;

	for(int ray=0; ray<nRays; ray++)
	{
		if(m_steps[ray] == -1) return;
		else if(m_steps[ray] >= STEPLIM) return;

		int n;
		n = Propagate(ray, r_max, theta_max, STEPLIM, outfile, write_step, write_rmax, write_rmin, write_cartesian);
		m_steps[ray] += n;

		if(outfile != 0)
			outfile->newline(2);
	}
}

template <typename T>
inline int Raytracer<T>::Propagate(int ray, const T rlim, const T thetalim, const int steplim, TextOutput* outfile, int write_step, T write_rmax, T write_rmin, bool write_cartesian )
{
	//
	// propagate the photon along its geodesic until limiting r or theta reached
	//
	int steps = 0;

	int rsign_count = COUNT_MIN;
	int thetasign_count = COUNT_MIN;

	T rhosq, delta, sigmasq, e2nu, e2psi, omega, grr, gthth, gphph;
	T rdotsq, thetadotsq;
	T energy;
	T len;

	T step;

	T x, y, z;

	// copy variables locally to simplify the code
	T a = spin;
	T t = m_t[ray];
	T r = m_r[ray];
	T theta = m_theta[ray];
	T phi = m_phi[ray];
	T pt = m_pt[ray];
	T pr = m_pr[ray];
	T ptheta = m_ptheta[ray];
	T pphi = m_pphi[ray];
	int rdot_sign = m_rdot_sign[ray];
	int thetadot_sign = m_thetadot_sign[ray];

	const T k = m_k[ray];
	const T h = m_h[ray];
	const T Q = m_Q[ray];

	bool write_started = false;

	// integrate geodesic equations until limit reached
	// if thetalim is positive, we go until theta exceeds it, if it is negative, we go until it is less than the abs value to allow tracing back to theta=0
	while( r < rlim  && ( (thetalim > 0 && theta < thetalim) || (thetalim < 0 && theta > abs(thetalim)) || thetalim == 0 )  &&  steps < steplim )
	//while( r < rlim  && theta < thetalim  &&  steps < steplim )
	{
		++steps;

		rhosq = r*r + (a*cos(theta))*(a*cos(theta));
		delta = r*r - 2*r + a*a;
		sigmasq = (r*r + spin*spin)*(r*r + spin*spin) - spin*spin*delta*sin(theta)*sin(theta);
		e2nu = rhosq * delta / sigmasq;
		e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
		omega = 2*spin*r / sigmasq;

		// tdot
		pt = (rhosq*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta))*k - 2*a*r*h;
		pt = pt / ( r*r * (1 + (a*cos(theta)/r)*(a*cos(theta)/r) - 2/r)*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta) );

		// phidot
		pphi = 2*a*r*sin(theta)*sin(theta)*k + (r*r + (a*cos(theta))*(a*cos(theta)) - 2*r)*h;
		pphi = pphi / ( (r*r + a*a)*(r*r + (a*cos(theta))*(a*cos(theta)) - 2*r)*sin(theta)*sin(theta) + 2*a*a*r*sin(theta)*sin(theta)*sin(theta)*sin(theta) );

		// thetadot
		thetadotsq = Q + (k*a*cos(theta) + h/tan(theta))*(k*a*cos(theta) - h/tan(theta));
		thetadotsq = thetadotsq / (rhosq*rhosq);

		if(thetadotsq < 0 && thetasign_count >= COUNT_MIN)
		{
			thetadot_sign *= -1;
			thetasign_count = 0;
			continue;
		}
		if (thetasign_count <= COUNT_MIN) thetasign_count++;

		// take the square roots and get the right signs
		ptheta = sqrt(abs(thetadotsq)) * thetadot_sign;

		// rdot
		rdotsq = k*pt - h*pphi - rhosq*ptheta*ptheta;
		rdotsq = rdotsq * delta/rhosq;

		if(rdotsq < 0 && rsign_count >= COUNT_MIN)
		{
			rdot_sign *= -1;
			rsign_count = 0;
			continue;
		}
		if (rsign_count <= COUNT_MIN) rsign_count++;

		pr = sqrt(abs(rdotsq)) * rdot_sign;

		step = abs( (r-(T)horizon)/pr ) / tolerance;
		// if the step is smaller in theta or phi (near 0/pi/2pi), use that instead
		if( step > abs( theta/ptheta ) / tolerance ) step = abs( theta/ptheta ) / tolerance;
		if( step > abs( phi/pphi ) / tolerance ) step = abs( phi/pphi ) / tolerance;
//		if( step > abs( (phi - M_PI)/pphi ) / tol ) step = abs( (phi - M_PI)/pphi ) / tol;
//		if( step > abs( (phi - 2*M_PI)/pphi ) / tol ) step = abs( (phi - 2*M_PI)/pphi ) / tol;
		// don't let the step be stupidly small
		if( step < MIN_STEP ) step = MIN_STEP;

		// calculate new position
		T dt += pt*step;
		T dr += pr*step;
		T dtheta += ptheta*step;
		T dphi += pphi*step;
		
		t += dt;
		r += dr;
		theta += dtheta;
		phi += dphi;

		grr = -rhosq/delta;
		gthth = -rhosq;
		gphph = -e2psi;
		
		len = grr*dr*dr + gthth*dtheta*dtheta + gphph*dphi*dphi;
		dens = 1;

		energy = 1./ray_redshift(V, false, false, r, theta, phi, k, h, Q, rdot_sign, thetadot_sign, m_emit[ray]);
		int ien = (logbin_en) ? static_cast<int>( log(energy/en0) / log(den) ) : static_cast<int>((energy - en0)/den);

		emis[m_ray][ien] += len * dens * pow(energy, 3);
		absorb[m_ray][ien] += len * dens;

		//aff += step;

		if(r <= horizon) break;

		if(outfile != 0 && (steps % write_step) == 0 )
		{
			if((write_rmax < 0 || r < write_rmax) && (write_rmin < 0 || r > write_rmin) )
			{
				write_started = true;
				if(write_cartesian)
				{
					Cartesian<T>(x, y, z, r, theta, phi, a);
					(*outfile) << t << x << y << z << endl;
				}
				else
				{
					(*outfile) << t << r << theta << phi << endl;
				}
			}
			else if(write_started)
			{
				break;
			}
		}
	}

	m_t[ray] = t;
	m_r[ray] = r;
	m_theta[ray] = theta;
	m_phi[ray] = phi;
	m_pt[ray] = pt;
	m_pr[ray] = pr;
	m_ptheta[ray] = ptheta;
	m_pphi[ray] = pphi;
	m_rdot_sign[ray] = rdot_sign;
	m_thetadot_sign[ray] = thetadot_sign;

	return steps;
}

template <typename T>
void Raytracer<T>::RedshiftStart( T V, bool reverse, bool projradius )
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
	// Arguments:
	//	V		T		Angular velocity of emitter (set to -1 for circular orbit at current r of ray)
	//	reverse	bool	Whether the ray is being propagated backwards in time (for image planes, default value = false)
	//
	cout << "Calculating initial energies" << endl;

	for(int ray=0; ray<nRays; ray++)
	{
		T p[4];

		const T a = (reverse) ? -1*spin : spin;

		// metric coefficients
		const T rhosq = m_r[ray]*m_r[ray] + (a*cos(m_theta[ray]))*(a*cos(m_theta[ray]));
		const T delta = m_r[ray]*m_r[ray] - 2*m_r[ray] + a*a;
		const T sigmasq = (m_r[ray]*m_r[ray] + a*a)*(m_r[ray]*m_r[ray] + a*a) - a*a*delta*sin(m_theta[ray])*sin(m_theta[ray]);

		const T e2nu = rhosq * delta / sigmasq;
		const T e2psi = sigmasq * sin(m_theta[ray])*sin(m_theta[ray]) / rhosq;
		const T omega = 2*a*m_r[ray] / sigmasq;

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

		// if V==-1, calculate orbital velocity for a geodesic circular orbit in equatorial plane
		if(V == -1 && projradius)
			V = 1 / (a + m_r[ray]*sin(m_theta[ray])*sqrt(m_r[ray]*sin(m_theta[ray])));	// project the radius parallel to the equatorial plane
		else if(V == -1)
			V = 1 / (a + m_r[ray]*sqrt(m_r[ray]));


		// timelike basis vector
		const T et[] = { (1/sqrt(e2nu))/sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu)
							, 0 , 0 ,
							(1/sqrt(e2nu))*V / sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu) };

		// photon momentum
		MomentumFromConsts<T>(p[0], p[1], p[2], p[3], m_k[ray], m_h[ray], m_Q[ray], m_rdot_sign[ray], m_thetadot_sign[ray], m_r[ray], m_theta[ray], m_phi[ray], spin);

		// if we're propagating backwards, reverse the direction of the photon momentum
		if(reverse) p[1] *= -1; p[2] *= -1; p[3] *= -1;

		// evaluate dot product to get energy
		m_emit[ray] = 0;
		for(int i=0; i<4; i++)
			for(int j=0; j<4; j++)
				m_emit[ray] += g[i*4 + j] * et[i]* p[j];
	}
}


template <typename T>
void Raytracer<T>::Redshift( T V, bool reverse, bool projradius )
{
	//
	// Calculates the redshift of the ray (emitted / received energy).
	// Receiver is orbitting the black hole rotation axis at angular velocity V, or if V = -1, the
	// angular velocity is calculated for a circular orbit at the current r co-ordinate of the ray.
	//
	// If ray is being propagated backwards in time (for image planes), the spin parameter is revered
	// (back to its true value as it will already have been reversed for the propagation) and the
	// spatial components of the ray's 4-momentum will be reversed so the ray is travelling in the
	// correct direction wrt the receiving material.
	//
	// RedshiftStart( ) needs to have been called before rays were traced. Use this function after
	// propagation.
	//
	// Arguments:
	//	V		T		Angular velocity of receiver (set to -1 for circular orbit at current r of ray)
	//	reverse	bool	Whether the ray is being propagated backwards in time (for image planes, default value = false)
	//	projradius bool	Whether to use r*sin(theta) as the radial co-ordinate in angular velocity calculation to use radius projected parallel to equatorial plane
	//
	cout << "Calculating ray redshifts..." << endl;

	for(int ray=0; ray<nRays; ray++)
	{
		T p[4];

		// metric coefficients
		const T rhosq = m_r[ray]*m_r[ray] + (spin*cos(m_theta[ray]))*(spin*cos(m_theta[ray]));
		const T delta = m_r[ray]*m_r[ray] - 2*m_r[ray] + spin*spin;
		const T sigmasq = (m_r[ray]*m_r[ray] + spin*spin)*(m_r[ray]*m_r[ray] + spin*spin) - spin*spin*delta*sin(m_theta[ray])*sin(m_theta[ray]);

		const T e2nu = rhosq * delta / sigmasq;
		const T e2psi = sigmasq * sin(m_theta[ray])*sin(m_theta[ray]) / rhosq;
		const T omega = 2*spin*m_r[ray] / sigmasq;

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

		// if V==-1, calculate orbital velocity for a geodesic circular orbit in equatorial lane
		if(V == -1 && projradius)
			V = 1 / (spin + m_r[ray]*sin(m_theta[ray])*sqrt(m_r[ray]*sin(m_theta[ray])));	// project the radius parallel to the equatorial plane
		else if(V == -1)
			V = 1 / (spin + m_r[ray]*sqrt(m_r[ray]));

		// timelike basis vector
		const T et[] = { (1/sqrt(e2nu))/sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu)
							, 0 , 0 ,
							(1/sqrt(e2nu))*V / sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu) };

		// photon momentum
		MomentumFromConsts<T>(p[0], p[1], p[2], p[3], m_k[ray], m_h[ray], m_Q[ray], m_rdot_sign[ray], m_thetadot_sign[ray], m_r[ray], m_theta[ray], m_phi[ray], spin);

		// if we're propagating backwards, reverse the direction of the photon momentum
		if(reverse) p[1] *= -1; p[2] *= -1; p[3] *= -1;

		// evaluate dot product to get energy
		T recv = 0;
		for(int i=0; i<4; i++)
			for(int j=0; j<4; j++)
				recv += g[i*4 + j] * et[i] * p[j];

		if(reverse)
			m_redshift[ray] = recv / m_emit[ray];
		else
			m_redshift[ray] = m_emit[ray] / recv;
	}
}



template <typename T>
void Raytracer<T>::RangePhi( T min, T max )
{
	//
	// Puts the azimuthal angle co-ordinate (phi) in the required range
	//
	// Arguments:
	//	min		T		Lower bound of range (default value = -1 * M_PI)
	//	max		T		Upper bound of range (default value = M_PI)
	//
	cout << "Putting phi co-ordinate into range [" << min << "," << max << "]" << endl;
	for( int ray=0; ray<nRays; ray++ )
	{
		// check phi isn't something horrible so we don't enter an infinite loop
		if( abs(m_phi[ray]) > 1000 || m_phi[ray] != m_phi[ray] || !(m_steps[ray]>0) ) continue;

		while( m_phi[ray] >= max ) m_phi[ray] -= 2*M_PI;
		while( m_phi[ray] < min ) m_phi[ray] += 2*M_PI;
	}
}


template <typename T>
inline void Raytracer<T>::CalculateConstants(int ray, T alpha, T beta, T V, T E)
{
	//
	// Compute the constants of motion for a ray emitted at polar angles alpha and beta in the frame
	// of a source at (m_t[ray],m_r[ray],m_theta[ray],m_phi[ray]) orbiting azimuthally at angular velocity V
	//
	const T rhosq = m_r[ray]*m_r[ray] + (spin*cos(m_theta[ray]))*(spin*cos(m_theta[ray]));
	const T delta = m_r[ray]*m_r[ray] - 2*m_r[ray] + spin*spin;
	const T sigmasq = (m_r[ray]*m_r[ray] + spin*spin)*(m_r[ray]*m_r[ray] + spin*spin) - spin*spin*delta*sin(m_theta[ray])*sin(m_theta[ray]);

	// metric coefficients
	const T e2nu = rhosq * delta / sigmasq;
	const T e2psi = sigmasq * sin(m_theta[ray])*sin(m_theta[ray]) / rhosq;
	const T omega = 2*spin*m_r[ray] / sigmasq;

	// tetrad basis vector components
	const T et0 = (1/sqrt(e2nu))/sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu);
	const T et3 = (1/sqrt(e2nu))*V / sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu);
	//
	const T e10 = (V - omega)*sqrt(e2psi/e2nu) / sqrt(e2nu - (V-omega)*(V-omega)*e2psi);
	const T e13 = (1/sqrt(e2nu*e2psi))*(e2nu + V*omega*e2psi - omega*omega*e2psi) / sqrt(e2nu - (V-omega)*(V-omega)*e2psi);
	//
	const T e22 = 1/sqrt(rhosq);
	//
	const T e31 = sqrt(delta/rhosq);

	// photon 4-momentum in source frame
	const T rdotprime[] = { E, E*sin(alpha)*sin(beta), E*sin(alpha)*cos(beta), E*cos(alpha) };

	const T tdot = rdotprime[0]*et0 + rdotprime[1]*e10;
	const T phidot = rdotprime[0]*et3 + rdotprime[1]*e13;
	// const T rdot = rdotprime[3]*e31;
	const T thetadot = rdotprime[2]*e22;

	// find the corresponding values of k, h and Q using the geodesic equations
	m_k[ray] = (1 - 2*m_r[ray]/rhosq)*tdot + (2*spin*m_r[ray]*sin(m_theta[ray])*sin(m_theta[ray])/rhosq)*phidot;

	m_h[ray] = phidot * ( (m_r[ray]*m_r[ray] + spin*spin)*(m_r[ray]*m_r[ray] + spin*spin*cos(m_theta[ray])*cos(m_theta[ray]) - 2*m_r[ray])*sin(m_theta[ray])*sin(m_theta[ray]) + 2*spin*spin*m_r[ray]*sin(m_theta[ray])*sin(m_theta[ray])*sin(m_theta[ray])*sin(m_theta[ray]) );
	m_h[ray] = m_h[ray] - 2*spin*m_r[ray]*m_k[ray]*sin(m_theta[ray])*sin(m_theta[ray]);
	m_h[ray] = m_h[ray] / ( m_r[ray]*m_r[ray] + spin*spin*cos(m_theta[ray])*cos(m_theta[ray]) - 2*m_r[ray] );

	m_Q[ray] = rhosq*rhosq*thetadot*thetadot - (spin*m_k[ray]*cos(m_theta[ray]) + m_h[ray]/tan(m_theta[ray]))*(spin*m_k[ray]*cos(m_theta[ray]) - m_h[ray]/tan(m_theta[ray]));
}

template <typename T>
inline void Raytracer<T>::CalculateConstantsFromP(int ray, T pt, T pr, T ptheta, T pphi)
{
	//
	// calculate constants of motion from the 4-momentum and location of a photon
	//
	const T a = spin;
	const T r = m_r[ray];
	const T theta = m_theta[ray];

	const T rhosq = r*r + (a*cos(theta))*(a*cos(theta));

	T k = (1 - 2*r/rhosq)*pr + (2*a*r*sin(theta)*sin(theta)/rhosq)*pphi;

	T h = pphi * ( (r*r + a*a)*(r*r + a*a*cos(theta)*cos(theta) - 2*r)*sin(theta)*sin(theta) + 2*a*a*r*sin(theta)*sin(theta)*sin(theta)*sin(theta) );
	h = h - 2*a*r*k*sin(theta)*sin(theta);
	h = h / ( r*r + a*a*cos(theta)*cos(theta) - 2*r );

	T Q = rhosq*rhosq*ptheta*ptheta - (a*k*cos(theta) + h/tan(theta))*(a*k*cos(theta) - h/tan(theta));
	
	m_k[ray] = k;
	m_h[ray] = h;
	m_Q[ray] = Q;
}


template<typename T>
void Raytracer<T>::CalculateMomentum( )
{
	//
	// Calculate photon momentum from constants of motion at a location
	//
	T a = spin;

	for(int ray = 0; ray < nRays; ray++)
	{
		const T t = m_t[ray];
		const T r = m_r[ray];
		const T theta = m_theta[ray];
		const T phi = m_phi[ray];
		const int rdot_sign = m_rdot_sign[ray];
		const int thetadot_sign = m_thetadot_sign[ray];

		const T k = m_k[ray];
		const T h = m_h[ray];
		const T Q = m_Q[ray];

		const T rhosq = r*r + (a*cos(theta))*(a*cos(theta));
		const T delta = r*r - 2*r + a*a;

		// tdot
		m_pt[ray] = (rhosq*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta))*k - 2*a*r*h;
		m_pt[ray] = m_pt[ray] / ( r*r * (1 + (a*cos(theta)/r)*(a*cos(theta)/r) - 2/r)*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta) );

		// phidot
		m_pphi[ray] = 2*a*r*sin(theta)*sin(theta)*k + (r*r + (a*cos(theta))*(a*cos(theta)) - 2*r)*h;
		m_pphi[ray] = m_pphi[ray] / ( (r*r + a*a)*(r*r + (a*cos(theta))*(a*cos(theta)) - 2*r)*sin(theta)*sin(theta) + 2*a*a*r*sin(theta)*sin(theta)*sin(theta)*sin(theta) );

		// thetadot
		T thetadotsq = Q + (k*a*cos(theta) + h/tan(theta))*(k*a*cos(theta) - h/tan(theta));
		thetadotsq = thetadotsq / (rhosq*rhosq);

		// take the square roots and get the right signs
		m_ptheta[ray] = sqrt(abs(thetadotsq)) * thetadot_sign;

		// rdot
		T rdotsq = k*m_pt[ray] - h*m_pphi[ray] - rhosq*m_ptheta[ray]*m_ptheta[ray];
		rdotsq = rdotsq * delta/rhosq;

		m_pr[ray] = sqrt(abs(rdotsq)) * rdot_sign;
	}
}

template class Raytracer<double>;
template class Raytracer<float>;
