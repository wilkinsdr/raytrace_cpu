/*
 * imagePlane.cpp
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#include "sourcetracer_imageplane.h"

template <typename T>
SourceTracer_ImagePlane<T>::SourceTracer_ImagePlane( T dist, T inc, T x0, T xmax, T dx, T y0, T ymax, T dy, T spin, T init_en0, T init_enmax, int init_Nen, bool init_logbin_en, T tol, T phi )
	: SourceTracer<T>( (((xmax - x0) / dx) + 1) * (((ymax - y0) / dy) + 1) , spin , init_en0, init_enmax, init_Nen, init_logbin_en, tol, true ),
	        Nx(((xmax - x0) / dx) + 1),
	        Ny(((ymax - y0) / dy) + 1),
	        D(dist),
	        incl(inc),
	        phi0(phi)
{
	m_plane_x = new T[Raytracer<T>::nRays];
	m_plane_y = new T[Raytracer<T>::nRays];

	cout << "Setting up image plane with (" << Nx << 'x' << Ny << ") rays" << endl;
	InitImagePlane( D, incl*M_PI/180, phi0, x0, xmax, dx, y0, ymax, dy);
}

template <typename T>
SourceTracer_ImagePlane<T>::~SourceTracer_ImagePlane()
{
	delete[] m_plane_x;
	delete[] m_plane_y;
}

template <typename T>
void SourceTracer_ImagePlane<T>::InitImagePlane( T D, T incl, T phi0,
									T x0, T xmax, T dx,
                                    T y0, T ymax, T dy)
{
	const int Nx = ((xmax - x0) / dx) + 1;
	const int Ny = ((ymax - y0) / dy) + 1;
	
	const double a = Raytracer<T>::spin;
	
	for(int i=0; i<Nx; i++)
	{
		const T x = x0 + i * dy;

		for (int j = 0; j < Ny; j++)
		{
			const int ix = i * Ny + j;

			const T y = y0 + j*dy;

			// initialise position of photon
			const T t = 0;
			const T r = sqrt(D*D + x*x + y*y);
			const T theta = acos( (D*cos(incl) + y*sin(incl)) / r );
			const T phi = phi0 + atan2( x, D*sin(incl) - y*cos(incl) );
			
			// and the momentum
			const T pr = D/r;
			const T ptheta = sin(acos(D/r))/r;
			const T pphi = x*sin(incl)/(x*x+(D*sin(incl)-y*cos(incl))*(D*sin(incl)-y*cos(incl)));

			// metric coefficients
			const T rhosq = r*r + (a*cos(theta))*(a*cos(theta));
			const T delta = r*r - 2*r + a*a;
			const T sigmasq = (r*r + a*a)*(r*r + a*a) - a*a*delta*sin(theta)*sin(theta);

			const T e2nu = rhosq * delta / sigmasq;
			const T e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
			const T omega = 2*a*r / sigmasq;

			const T g00 = e2nu - omega*omega*e2psi;
			const T g03 = omega*e2psi;
			const T g11 = -rhosq/delta;
			const T g22 = -rhosq;
			const T g33 = -e2psi;

			// solve quadratic in pt to make this a null vector
			const T A = g00;
			const T B = 2*g03*pphi;
			const T C = g11*pr*pr + g22*ptheta*ptheta + g33*pphi*pphi;

			// take the appropriate (positive) root
			T pt = ( -B + sqrt(B*B - 4*A*C) ) / (2*A);
			if(pt < 0) pt = ( -B - sqrt(B*B - 4*A*C) ) / (2*A);

			Raytracer<T>::m_t[ix] = t;
			Raytracer<T>::m_r[ix] = r;
			Raytracer<T>::m_theta[ix] = theta;
			Raytracer<T>::m_phi[ix] = phi;
			Raytracer<T>::m_pt[ix] = pt;
			Raytracer<T>::m_pr[ix] = pr;
			Raytracer<T>::m_ptheta[ix] = ptheta;
			Raytracer<T>::m_pphi[ix] = pphi;

			// calculate constants of motion
			Raytracer<T>::CalculateConstantsFromP(ix, pt, pr, ptheta, pphi);
			Raytracer<T>::m_rdot_sign[ix] = -1;
			Raytracer<T>::m_thetadot_sign[ix] = 1;

			Raytracer<T>::m_k[ix] = 1;

			const T b = sqrt(x*x + y*y);
			T beta = asin(y/b);
			if(x < 0) beta=M_PI-beta;

			T h = -1.*b*sin(incl)*cos(beta);
			T ltheta = b*sin(beta);
			T Q = (ltheta*ltheta) - (a*cos(theta))*(a*cos(theta))+((h/tan(theta)))*((h/tan(theta)));

			Raytracer<T>::m_h[ix] = h;
			Raytracer<T>::m_Q[ix] = Q;

			Raytracer<T>::m_thetadot_sign[ix] = (ltheta>=0) ? 1 : -1;
			
			Raytracer<T>::m_steps[ix] = 0;
		}
	}
}

template <typename T>
void SourceTracer_ImagePlane<T>::RedshiftStart( )
{
	//
	// Call the RedshiftStart function of the base class using the source's angular velocity
	//
	Raytracer<T>::RedshiftStart( 0, true );
}

template <typename T>
void SourceTracer_ImagePlane<T>::Redshift( bool projradius )
{
	//
	// Call the RedshiftStart function of the base class using the angular velocity for a circular orbit at the ray's end point
	// for rays incident on the accretion disc
	//
	Raytracer<T>::Redshift( -1, true, projradius );
}

template class SourceTracer_ImagePlane<double>;
template class SourceTracer_ImagePlane<float>;
