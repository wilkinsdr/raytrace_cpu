/*
 * pointsource_vel.cpp
 *
 *  Created on: 13 Jul 2017
 *      Author: drw
 */

#include "pointsource_vel.h"

template <typename T>
PointSourceVel<T>::PointSourceVel( T* pos, T* V, T spin, T tol, T dcosalpha, T dbeta, T cosalpha0, T cosalphamax, T beta0, T betamax, T E )
	: Raytracer<T>( (((cosalphamax - cosalpha0) / dcosalpha) + 1) * (((betamax - beta0) / dbeta) + 1) , spin , tol ),
	  energy(E)
{
	n_cosalpha = ((cosalphamax - cosalpha0) / dcosalpha) + 1;
	n_beta = ((betamax - beta0) / dbeta) + 1;

	m_cosalpha = new T[Raytracer<T>::nRays];
	m_beta = new T[Raytracer<T>::nRays];

	cout << "Setting up point source with " << Raytracer<T>::nRays << " rays" << endl;
    SolveTetrad( pos, V, spin );
	InitPointSourceVel( pos, dcosalpha, dbeta, cosalpha0, cosalphamax, beta0, betamax );
}

template <typename T>
PointSourceVel<T>::~PointSourceVel()
{
	delete[] m_cosalpha;
	delete[] m_beta;
}

template <typename T>
void PointSourceVel<T>::InitPointSourceVel( T* pos, T dcosalpha, T dbeta, T cosalpha0, T cosalphamax, T beta0, T betamax )
{
	for(int i=0; i<n_cosalpha; i++)
		for(int j=0; j<n_beta; j++)
		{
			int ix = i*n_beta + j;

            const T cosalpha = cosalpha0 + i*dcosalpha;
            const T beta = beta0 + j*dbeta;

			m_cosalpha[ix] = cosalpha;
			m_beta[ix] = beta;

			if(m_cosalpha[ix] >= cosalphamax || m_beta[ix] >= betamax)
			{
				Raytracer<T>::m_steps[ix] = -1;
				return;
			}

			const T alpha = acos(m_cosalpha[ix]);

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

            const T r = pos[1];
            const T theta = pos[2];
            const T phi = pos[3];
            const T a = Raytracer<T>::spin;

            //
            // Compute the constants of motion for a ray emitted at polar angles alpha and beta in the frame
            // of a source at (m_t[ray],r,theta,m_phi[ray]) using the tetrad calculated for the
            // source velocity
            //
            const T rhosq = r*r + (a*cos(theta))*(a*cos(theta));
            const T delta = r*r - 2*r + a*a;
            const T sigmasq = (r*r + a*a)*(r*r + a*a) - a*a*delta*sin(theta)*sin(theta);

            // metric coefficients
            const T e2nu = rhosq * delta / sigmasq;
            const T e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
            const T omega = 2*a*r / sigmasq;

            // photon 4-momentum in source frame
            const T rdotprime[] = { energy, energy*sin(alpha)*cos(beta), energy*sin(alpha)*sin(beta), energy*cos(alpha) };

            const T tdot = rdotprime[0]*m_tetrad[0][0] + rdotprime[1]*m_tetrad[1][0] + rdotprime[2]*m_tetrad[2][0] + rdotprime[3]*m_tetrad[3][0];
            const T rdot = rdotprime[0]*m_tetrad[0][1] + rdotprime[1]*m_tetrad[1][1] + rdotprime[2]*m_tetrad[2][1] + rdotprime[3]*m_tetrad[3][1];
            const T thetadot = rdotprime[0]*m_tetrad[0][2] + rdotprime[1]*m_tetrad[1][2] + rdotprime[2]*m_tetrad[2][2] + rdotprime[3]*m_tetrad[3][2];
            const T phidot = rdotprime[0]*m_tetrad[0][3] + rdotprime[1]*m_tetrad[1][3] + rdotprime[2]*m_tetrad[2][3] + rdotprime[3]*m_tetrad[3][3];

            // find the corresponding values of k, h and Q using the geodesic equations
            double k = (1 - 2*r/rhosq)*tdot + (2*a*r*sin(theta)*sin(theta)/rhosq)*phidot;

            double h = phidot * ( (r*r + a*a)*(r*r + a*a*cos(theta)*cos(theta) - 2*r)*sin(theta)*sin(theta) + 2*a*a*r*sin(theta)*sin(theta)*sin(theta)*sin(theta) );
            h = h - 2*a*r*k*sin(theta)*sin(theta);
            h = h / ( r*r + a*a*cos(theta)*cos(theta) - 2*r );

            double Q = rhosq*rhosq*thetadot*thetadot - (a*k*cos(theta) + h/tan(theta))*(a*k*cos(theta) - h/tan(theta));

            Raytracer<T>::m_k[ix] = k;
            Raytracer<T>::m_h[ix] = h;
            Raytracer<T>::m_Q[ix] = Q;

            Raytracer<T>::m_rdot_sign[ix] = (rdot > 0) ? 1 : -1;
            Raytracer<T>::m_thetadot_sign[ix] = (thetadot > 0) ? 1 : -1;

		}
}

template <typename T>
void PointSourceVel<T>::SolveTetrad(T* pos, T* vel, T spin)
{
	T g[4][4];
    T v[4][4];
    T e[4][4];

    const T t = pos[0];
	const T r = pos[1];
	const T theta = pos[2];
	const T phi = pos[3];

    const double rhosq = r*r + (spin*cos(theta))*(spin*cos(theta));
    const double delta = r*r - 2*r + spin*spin;
    const double sigmasq = (r*r + spin*spin)*(r*r + spin*spin) - spin*spin*delta*sin(theta)*sin(theta);

    const double e2nu = rhosq * delta / sigmasq;
    const double e2psi = sigmasq * sin(theta)*sin(theta) / rhosq;
    const double omega = 2*spin*r / sigmasq;

    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
        {
            g[i][j] = 0;
            v[i][j] = 0;
            e[i][j] = 0;
        }

    // metric co-efficients
    g[0][0] = e2nu - omega*omega*e2psi;
    g[0][3] = omega*e2psi;
    g[3][0] = g[0][3];
    g[1][1] = -rhosq/delta;
    g[2][2] = -rhosq;
    g[3][3] = -e2psi;

    // initial guess of the basis vectors
    // v0 is the 4-velocity
    v[0][0] = vel[0];
    v[0][1] = vel[1];
    v[0][2] = vel[2];
    v[0][3] = vel[3];
    // 1, 2 and 3 are initially in the phi, theta and r directions
    v[1][1] = 1; // put the radial one first to prioritise getting one in approximately this direction
    v[2][2] = 1;
    v[3][3] = 1;

    cout << "4-velocity:" << endl << setw(15) << v[0][0] << setw(15) << v[0][1] << setw(15) << v[0][2] << setw(15) << v[0][3] << endl;

    // u0 = v0
    for(int a=0; a<4; a++)
        e[0][a] = v[0][a];

    // Gram-Schmidt orthogonalisation
    for(int i=1; i<4; i++)
    {
        for(int a=0; a<4; a++)
        {
            e[i][a] = v[i][a];
        }

        for(int j=0; j<i; j++)
        {
            double dotprod = 0;
            double enorm = 0;

            // evaluate the dot products
            for(int a=0; a<4; a++)
                for(int b=0; b<4; b++)
                {
                    dotprod += g[a][b]*v[i][a]*e[j][b];
                    enorm += g[a][b]*e[j][a]*e[j][b];
                }
            // subtract the projections of the existing u's from this one
            for(int a=0; a<4; a++)
            {
                e[i][a] -= (dotprod/enorm) * e[j][a];
            }
        }
    }

    cout << endl << "Orthogonalised vectors:" << endl;
    cout << setw(18) << 0 << setw(15) << 1 << setw(15) << 2 << setw(15) << 3 << endl;
    for(int i=0; i<4; i++)
    {
        cout << i << ": ";
        for(int a=0; a<4; a++)
            cout << setw(15) << e[i][a];
        cout << endl;
    }

    // normalise the vectors
    for(int i=0; i<4; i++)
    {
        double enorm = 0;

        // evaluate the dot product
        for(int a=0; a<4; a++)
            for(int b=0; b<4; b++)
                enorm += g[a][b]*e[i][a]*e[i][b];

        for(int a=0; a<4; a++)
            e[i][a] /= sqrt(abs(enorm));
    }

    // and store them in the object to use later
    // we need to do some swapping to make a right-handed set
    for(int i=0; i<4; i++) {
        int iout;
        switch(i) {
            case 0:
                iout = 0;
                break;
            case 1:
                iout = 3;
                break;
            case 2:
                iout = 1;
                break;
            case 3:
                iout = 2;
                break;
            default:
                iout = i;
        }
        for (int j = 0; j < 4; j++) {
            m_tetrad[iout][j] = e[i][j];
        }
    }

    // make sure the vectors are also pointing the correct way
    // for a right-handed basis
    if(m_tetrad[1][2] < 0)
        for(int i=0; i<4; i++) m_tetrad[1][i] *= -1;
    if(m_tetrad[2][3] < 0)
        for(int i=0; i<4; i++) m_tetrad[2][i] *= -1;
    if(m_tetrad[3][1] < 0)
        for(int i=0; i<4; i++) m_tetrad[3][i] *= -1;

    cout << endl << "Tetrad:" << endl;
    cout << setw(18) << 0 << setw(15) << 1 << setw(15) << 2 << setw(15) << 3 << endl;
    for(int i=0; i<4; i++)
    {
        cout << i << ": ";
        for(int a=0; a<4; a++)
            cout << setw(15) << m_tetrad[i][a];
        cout << endl;
    }
}

template <typename T>
void PointSourceVel<T>::RedshiftStart( )
{
	//
	// Call the RedshiftStart function of the base class using the source's angular velocity
	//
	Raytracer<T>::RedshiftStart( velocity );
}

template <typename T>
void PointSourceVel<T>::Redshift( T V )
{
	//
	// Call the RedshiftStart function of the base class using the angular velocity for a circular orbit at the ray's end point
	// for rays incident on the accretion disc
	//
	Raytracer<T>::Redshift( V );
}

template class PointSourceVel<double>;
template class PointSourceVel<float>;
