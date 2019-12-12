/*
 * raytracer.cu
 *
 *  Created on: 4 Sep 2013
 *      Author: drw
 */

#include "source_tracer.h"

template <typename T>
SourceTracer<T>::SourceTracer( int num_rays, float spin_par, T init_en0, T init_enmax, int init_Nen, bool init_logbin_en, float toler, bool init_reverse )
	: Raytracer<T>(num_rays, spin_par, toler), en0(init_en0), enmax(init_enmax), Nen(init_Nen), logbin_en(init_logbin_en), reverse(init_reverse)
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

	den = (logbin_en) ? exp(log(enmax/en0)/(Nen-1)) : (enmax - en0)/(Nen - 1);
	energy = new T[Nen];
	for(int ien=0; ien<Nen; ien++)
		energy[ien] = (logbin_en) ? en0 * pow(den, ien) : en0 + den*ien;
}

template <typename T>
SourceTracer<T>::~SourceTracer( )
{
	for(int i=0; i<Nen; i++)
	{
		delete[] emis[i];
		delete[] absorb[i];
	}
	delete[] emis;
	delete[] absorb;
}

template <typename T>
void SourceTracer<T>::run_source_trace( T r_max, T theta_max, TextOutput* outfile, int write_step, T write_rmax, T write_rmin, bool write_cartesian )
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
	cout << "Running source raytracer..." << endl;

	for(int ray=0; ray<Raytracer<T>::nRays; ray++)
	{
		if(ray % 100 == 0) cout << "Ray " << ray << '/' << Raytracer<T>::nRays << endl;
		if(Raytracer<T>::m_steps[ray] == -1) return;
		else if(Raytracer<T>::m_steps[ray] >= STEPLIM) return;

		int n;
		n = propagate_source(ray, r_max, theta_max, STEPLIM, outfile, write_step, write_rmax, write_rmin, write_cartesian);
		Raytracer<T>::m_steps[ray] += n;

		if(outfile != 0)
			outfile->newline(2);
	}
}

template <typename T>
inline int SourceTracer<T>::propagate_source(int ray, const T rlim, const T thetalim, const int steplim, TextOutput* outfile, int write_step, T write_rmax, T write_rmin, bool write_cartesian )
{
	//
	// propagate the photon along its geodesic until limiting r or theta reached
	//
	int steps = 0;

	int rsign_count = COUNT_MIN;
	int thetasign_count = COUNT_MIN;
	
	T rdotsq, thetadotsq;

	T step;

	T x, y, z;

	// copy variables locally to simplify the code
	T a = Raytracer<T>::spin;
	T t = Raytracer<T>::m_t[ray];
	T r = Raytracer<T>::m_r[ray];
	T theta = Raytracer<T>::m_theta[ray];
	T phi = Raytracer<T>::m_phi[ray];
	T pt = Raytracer<T>::m_pt[ray];
	T pr = Raytracer<T>::m_pr[ray];
	T ptheta = Raytracer<T>::m_ptheta[ray];
	T pphi = Raytracer<T>::m_pphi[ray];
	int rdot_sign = Raytracer<T>::m_rdot_sign[ray];
	int thetadot_sign = Raytracer<T>::m_thetadot_sign[ray];

	const T k = Raytracer<T>::m_k[ray];
	const T h = Raytracer<T>::m_h[ray];
	const T Q = Raytracer<T>::m_Q[ray];

	bool write_started = false;

	// integrate geodesic equations until limit reached
	// if thetalim is positive, we go until theta exceeds it, if it is negative, we go until it is less than the abs value to allow tracing back to theta=0
	while( r < rlim  && ( (thetalim > 0 && theta < thetalim) || (thetalim < 0 && theta > abs(thetalim)) || thetalim == 0 )  &&  steps < steplim )
	//while( r < rlim  && theta < thetalim  &&  steps < steplim )
	{
		++steps;

		const T rhosq = r*r + (a*cos(theta))*(a*cos(theta));
		const T delta = r*r - 2*r + a*a;
		const T sigmasq = (r*r + a*a)*(r*r + a*a) - a*a*delta*sin(theta)*sin(theta);
		const T e2nu = rhosq * delta / sigmasq;
		const T e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
		const T omega = 2*a*r / sigmasq;
		const T grr = -rhosq/delta;
		const T gthth = -rhosq;
		const T gphph = -e2psi;

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

		step = abs( (r-(T)Raytracer<T>::horizon)/pr ) / Raytracer<T>::tolerance;
		// if the step is smaller in theta or phi (near 0/pi/2pi), use that instead
		if( step > abs( theta/ptheta ) / Raytracer<T>::tolerance ) step = abs( theta/ptheta ) / Raytracer<T>::tolerance;
		if( step > abs( phi/pphi ) / Raytracer<T>::tolerance ) step = abs( phi/pphi ) / Raytracer<T>::tolerance;
//		if( step > abs( (phi - M_PI)/pphi ) / tol ) step = abs( (phi - M_PI)/pphi ) / tol;
//		if( step > abs( (phi - 2*M_PI)/pphi ) / tol ) step = abs( (phi - 2*M_PI)/pphi ) / tol;
		// don't let the step be stupidly small
		if( step < MIN_STEP ) step = MIN_STEP;

		// calculate new position
		T dt = pt*step;
		T dr = pr*step;
		T dtheta = ptheta*step;
		T dphi = pphi*step;
		
		t += dt;
		r += dr;
		theta += dtheta;
		phi += dphi;

		if ( (r*r*sin(theta)*sin(theta)*cos(phi)*cos(phi)/(source_size_xy*source_size_xy) + r*r*sin(theta)*sin(theta)*sin(phi)*sin(phi)/(source_size_xy*source_size_xy) + r*r*cos(theta)*cos(theta)/(source_size_z*source_size_z)) < 1 )
		{
			const T len = -1*(grr * dr * dr + gthth * dtheta * dtheta + gphph * dphi * dphi);
			const T dens = 1;

			const T energy = 1. / Raytracer<T>::ray_redshift(source_vel, reverse, false, r, theta, phi, k, h, Q, rdot_sign, thetadot_sign, Raytracer<T>::m_emit[ray], source_motion);
			int ien = (logbin_en) ? static_cast<int>( log(energy / en0) / log(den)) : static_cast<int>((energy - en0) / den);

			if(ien >= 0 && ien < Nen)
			{
				emis[ray][ien] += (1./(r*r)) * len * dens * pow(energy, 3);
				absorb[ray][ien] += len * dens;
			}
		}

		if(r <= Raytracer<T>::horizon) break;

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

	Raytracer<T>::m_t[ray] = t;
	Raytracer<T>::m_r[ray] = r;
	Raytracer<T>::m_theta[ray] = theta;
	Raytracer<T>::m_phi[ray] = phi;
	Raytracer<T>::m_pt[ray] = pt;
	Raytracer<T>::m_pr[ray] = pr;
	Raytracer<T>::m_ptheta[ray] = ptheta;
	Raytracer<T>::m_pphi[ray] = pphi;
	Raytracer<T>::m_rdot_sign[ray] = rdot_sign;
	Raytracer<T>::m_thetadot_sign[ray] = thetadot_sign;

	return steps;
}

template class SourceTracer<double>;
template class SourceTracer<float>;
