/*
 * HealpixPointSource.cpp
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#include "healpix_pointsource.h"

template <typename T>
void calc_consts_from_vector(T& k, T& h, T& Q, int& rdot_sign, int& thetadot_sign, T E, T* vec, T& t, T& r, T& theta, T& phi, T V, T a, int motion = 0, int basis = 0)
{
    //
    // Compute the constants of motion for a ray emitted at polar angles alpha and beta in the frame
    // of a source at (t,r,theta,phi) orbiting azimuthally at angular velocity V
    //
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

    T tdot, rdot, thetadot, phidot;

    // photon 4-momentum in source frame
    const T rdotprime[] = { E, E*vec[0], E*vec[1], E*vec[2] };

    if(motion == 0)
    {
        // orbital motion
        if(basis == 0)
        {
            // tetrad basis vector components
            const T et0 = (1/sqrt(e2nu))/sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu);
            const T et3 = (1/sqrt(e2nu))*V / sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu);
            //
            const T e10 = (V - omega)*sqrt(e2psi/e2nu) / sqrt(e2nu - (V-omega)*(V-omega)*e2psi);
            const T e13 = (1/sqrt(e2nu*e2psi))*(e2nu + V*omega*e2psi - omega*omega*e2psi) / sqrt(e2nu - (V-omega)*(V-omega)*e2psi);
            //
            const T e22 = -1/sqrt(rhosq);
            //
            const T e31 = sqrt(delta/rhosq);

            tdot = rdotprime[0]*et0 + rdotprime[1]*e10;
            phidot = rdotprime[0]*et3 + rdotprime[1]*e13;
            rdot = rdotprime[3]*e31;
            thetadot = rdotprime[2]*e22;
        }
        else if(basis == 1)
        {
            // tetrad basis vector components
            const T et0 = (1/sqrt(e2nu))/sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu);
            const T et3 = (1/sqrt(e2nu))*V / sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu);
            //
            const T e10 = (V - omega)*sqrt(e2psi/e2nu) / sqrt(e2nu - (V-omega)*(V-omega)*e2psi);
            const T e13 = (1/sqrt(e2nu*e2psi))*(e2nu + V*omega*e2psi - omega*omega*e2psi) / sqrt(e2nu - (V-omega)*(V-omega)*e2psi);
            //
            const T e32 = -1/sqrt(rhosq);
            //
            const T e21 = -1 * sqrt(delta/rhosq);

            tdot = rdotprime[0]*et0 + rdotprime[1]*e10;
            phidot = rdotprime[0]*et3 + rdotprime[1]*e13;
            rdot = rdotprime[2]*e21;
            thetadot = rdotprime[3]*e32;
        }
    }
    else if(motion == 1)
    {
        // radial motion

        // tetrad basis vector components (they're right-handed!)
        const T et0 = 1 / sqrt(g00 + g11 * V * V);
        const T et1 = V * et0;

        const T e12 = 1 / sqrt(rhosq);

        const T e23 = sqrt(g00 / (g03 * g03 - g00 * g33));
        const T e20 = -(g03 / g00) * e23;

        const T e30 = sqrt(-g11 / g00) * V / sqrt(g00 + g11 * V * V);
        const T e31 = sqrt(-g00 / g11) / sqrt(g00 + g11 * V * V);

        tdot = rdotprime[0] * et0 + rdotprime[1] * e20 + rdotprime[3] * e30;
        phidot = rdotprime[1] * e23;
        rdot = rdotprime[0] * et1 + rdotprime[3] * e31;
        thetadot = rdotprime[2] * e12;
    }
    // find the corresponding values of k, h and Q using the geodesic equations
    k = (1 - 2*r/rhosq)*tdot + (2*a*r*sin(theta)*sin(theta)/rhosq)*phidot;

    h = phidot * ( (r*r + a*a)*(r*r + a*a*cos(theta)*cos(theta) - 2*r)*sin(theta)*sin(theta) + 2*a*a*r*sin(theta)*sin(theta)*sin(theta)*sin(theta) );
    h = h - 2*a*r*k*sin(theta)*sin(theta);
    h = h / ( r*r + a*a*cos(theta)*cos(theta) - 2*r );

    Q = rhosq*rhosq*thetadot*thetadot - (a*k*cos(theta) + h/tan(theta))*(a*k*cos(theta) - h/tan(theta));

    rdot_sign = (rdot > 0) ? 1 : -1;
    thetadot_sign = (thetadot > 0) ? 1 : -1;
}

template <typename T>
HealpixPointSource<T>::HealpixPointSource( T* pos, T V, T spin, int order, int motion, int basis, T tol)
	: Raytracer<T>( 5*12*(1<<order)*(1<<order) , spin , tol ),
	  velocity(V)
{
    nside = 1<<order;
    npix = 12*nside*nside;

	cout << "Setting up healpix point source with " << npix << " pixels (" << Raytracer<T>::nRays << " rays)" << endl;
	init_healpix_pointsource( pos, order, motion, basis );
}

template <typename T>
void HealpixPointSource<T>::init_healpix_pointsource(T* pos, int order, int motion, int basis )
{
    for(int pix=0; pix<npix; pix++)
    {
	    T corners[4*3];
	    T center[3];

        get_pixel_vectors<T>(order, pix, corners, center);

        // the corner rays
        for(int i=0; i<4; i++)
        {
            const int ix = 5*pix + i;

            // initialise position of photon
            Raytracer<T>::m_t[ix] = pos[0];
            Raytracer<T>::m_r[ix] = pos[1];
            Raytracer<T>::m_theta[ix] = pos[2];
            Raytracer<T>::m_phi[ix] = pos[3];

            Raytracer<T>::m_pt[ix] = 0;
            Raytracer<T>::m_pr[ix] = 0;
            Raytracer<T>::m_ptheta[ix] = 0;
            Raytracer<T>::m_pphi[ix] = 0;
            Raytracer<T>::m_steps[ix] = 0;

            // calculate constants of motion
            calc_consts_from_vector<T>(Raytracer<T>::m_k[ix], Raytracer<T>::m_h[ix], Raytracer<T>::m_Q[ix], Raytracer<T>::m_rdot_sign[ix], Raytracer<T>::m_thetadot_sign[ix], 1., &corners[3*i], Raytracer<T>::m_t[ix], Raytracer<T>::m_r[ix], Raytracer<T>::m_theta[ix], Raytracer<T>::m_phi[ix], velocity, Raytracer<T>::spin, motion, basis);
        }

        // and the centre ray
        Raytracer<T>::m_t[5*pix + 4] = pos[0];
        Raytracer<T>::m_r[5*pix + 4] = pos[1];
        Raytracer<T>::m_theta[5*pix + 4] = pos[2];
        Raytracer<T>::m_phi[5*pix + 4] = pos[3];

        Raytracer<T>::m_pt[5*pix + 4] = 0;
        Raytracer<T>::m_pr[5*pix + 4] = 0;
        Raytracer<T>::m_ptheta[5*pix + 4] = 0;
        Raytracer<T>::m_pphi[5*pix + 4] = 0;

        Raytracer<T>::m_steps[5*pix + 4] = 0;

        // calculate constants of motion
        calc_consts_from_vector<T>(Raytracer<T>::m_k[5*pix + 4], Raytracer<T>::m_h[5*pix + 4], Raytracer<T>::m_Q[5*pix + 4], Raytracer<T>::m_rdot_sign[5*pix + 4], Raytracer<T>::m_thetadot_sign[5*pix + 4], 1., center, Raytracer<T>::m_t[5*pix + 4], Raytracer<T>::m_r[5*pix + 4], Raytracer<T>::m_theta[5*pix + 4], Raytracer<T>::m_phi[5*pix + 4], velocity, Raytracer<T>::spin, motion, basis);

    }

}

template <typename T>
void HealpixPointSource<T>::redshift_start(bool reverse)
{
	//
	// Call the redshift_start function of the base class using the source's angular velocity
	//
    Raytracer<T>::redshift_start(velocity, reverse);
}

template <typename T>
void HealpixPointSource<T>::redshift(T V, bool reverse)
{
	//
	// Call the redshift_start function of the base class using the angular velocity for a circular orbit at the ray's end point
	// for rays incident on the accretion disc
	//
    Raytracer<T>::redshift(V, reverse);
}

template class HealpixPointSource<double>;
template class HealpixPointSource<float>;
