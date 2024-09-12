/*
 * imagePlane.cpp
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#include "mapper_imageplane.h"


Mapper_ImagePlane::Mapper_ImagePlane( double dist, double inc, double x0, double xmax, double dx, double y0, double ymax, double dy, double spin, double init_r0, double init_rmax, int init_Nr, int init_Ntheta, int init_Nphi, bool init_logbin_r, double init_thetamax, double tol, double phi )
	: Mapper( (((xmax - x0) / dx) + 1) * (((ymax - y0) / dy) + 1) , spin, init_r0, init_rmax, init_Nr, init_Ntheta, init_Nphi, init_logbin_r, init_thetamax, tol, true ),
	        Nx(((xmax - x0) / dx) + 1),
	        Ny(((ymax - y0) / dy) + 1),
	        D(dist),
	        incl(inc),
	        phi0(phi)
{
	m_plane_x = new double[Raytracer::nRays];
	m_plane_y = new double[Raytracer::nRays];

	cout << "Setting up image plane with (" << Nx << 'x' << Ny << ") rays" << endl;
	init_image_plane( D, incl*M_PI/180, phi0, x0, xmax, dx, y0, ymax, dy);
}


Mapper_ImagePlane::~Mapper_ImagePlane()
{
	delete[] m_plane_x;
	delete[] m_plane_y;
}


void Mapper_ImagePlane::init_image_plane( double D, double incl, double phi0,
									double x0, double xmax, double dx,
                                    double y0, double ymax, double dy)
{
	const int Nx = ((xmax - x0) / dx) + 1;
	const int Ny = ((ymax - y0) / dy) + 1;
	
	const double a = Raytracer::spin;
	
	for(int i=0; i<Nx; i++)
	{
		const double x = x0 + i * dy;

		for (int j = 0; j < Ny; j++)
		{
			const int ix = i * Ny + j;

			const double y = y0 + j*dy;

			// initialise position of photon
			const double t = 0;
			const double r = sqrt(D*D + x*x + y*y);
			const double theta = acos( (D*cos(incl) + y*sin(incl)) / r );
			const double phi = phi0 + atan2( x, D*sin(incl) - y*cos(incl) );
			
			// and the momentum
			const double pr = D/r;
			const double ptheta = sin(acos(D/r))/r;
			const double pphi = x*sin(incl)/(x*x+(D*sin(incl)-y*cos(incl))*(D*sin(incl)-y*cos(incl)));

			// metric coefficients
			const double rhosq = r*r + (a*cos(theta))*(a*cos(theta));
			const double delta = r*r - 2*r + a*a;
			const double sigmasq = (r*r + a*a)*(r*r + a*a) - a*a*delta*sin(theta)*sin(theta);

			const double e2nu = rhosq * delta / sigmasq;
			const double e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
			const double omega = 2*a*r / sigmasq;

			const double g00 = e2nu - omega*omega*e2psi;
			const double g03 = omega*e2psi;
			const double g11 = -rhosq/delta;
			const double g22 = -rhosq;
			const double g33 = -e2psi;

			// solve quadratic in pt to make this a null vector
			const double A = g00;
			const double B = 2*g03*pphi;
			const double C = g11*pr*pr + g22*ptheta*ptheta + g33*pphi*pphi;

			// take the appropriate (positive) root
			double pt = ( -B + sqrt(B*B - 4*A*C) ) / (2*A);
			if(pt < 0) pt = ( -B - sqrt(B*B - 4*A*C) ) / (2*A);

			Raytracer::m_t[ix] = t;
			Raytracer::m_r[ix] = r;
			Raytracer::m_theta[ix] = theta;
			Raytracer::m_phi[ix] = phi;
			Raytracer::m_pt[ix] = pt;
			Raytracer::m_pr[ix] = pr;
			Raytracer::m_ptheta[ix] = ptheta;
			Raytracer::m_pphi[ix] = pphi;

			// calculate constants of motion
			Raytracer::CalculateConstantsFromP(ix, pt, pr, ptheta, pphi);
			Raytracer::m_rdot_sign[ix] = -1;
			Raytracer::m_thetadot_sign[ix] = 1;

			Raytracer::m_k[ix] = 1;

			const double b = sqrt(x*x + y*y);
			double beta = asin(y/b);
			if(x < 0) beta=M_PI-beta;

			double h = -1.*b*sin(incl)*cos(beta);
			double ltheta = b*sin(beta);
			double Q = (ltheta*ltheta) - (a*cos(theta))*(a*cos(theta))+((h/tan(theta)))*((h/tan(theta)));

			Raytracer::m_h[ix] = h;
			Raytracer::m_Q[ix] = Q;

			Raytracer::m_thetadot_sign[ix] = (ltheta>=0) ? 1 : -1;
			
			Raytracer::m_steps[ix] = 0;
		}
	}
}


void Mapper_ImagePlane::redshift_start( )
{
	//
	// Call the redshift_start function of the base class using the source's angular velocity
	//
    Raytracer::redshift_start(0, true);
}


void Mapper_ImagePlane::redshift( bool projradius )
{
	//
	// Call the redshift_start function of the base class using the angular velocity for a circular orbit at the ray's end point
	// for rays incident on the accretion disc
	//
    Raytracer::redshift(-1, true, projradius);
}

