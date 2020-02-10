/*
 * kerr.h
 *
 *  Created on: 15 Oct 2013
 *      Author: drw
 */

#ifndef KERR_H_
#define KERR_H_

#include <cmath>

template<typename T>
T kerr_horizon(T a, int sign = 1)
{
	//
	// Calculate event horizon in Kerr geometry
	//
	return 1 + sign*sqrt( (1-a)*(1+a) );
}

template<typename T>
T kerr_isco(T a, int sign)
{
	//
	// Calculate radius of innermost stable circular orbit in Kerr geometry
	// for prograde (sign = +1) and retrograde (sign = -1) orbits
	//
	const float A = 1. + pow(1.-a*a , 1./3.) * ( pow(1.+a , 1./3.) + pow(1.-a , 1./3.) );
	const float B = sqrt(3.*a*a + A*A);
	return 3 + B - sign*sqrt( (3-A) * (3 + A + 2*B) );
}

template<typename T>
T disc_velocity(T r, T a, int sign)
{
	return 1 / (a + sign*pow(r, 3./2.));
}

template<typename T>
void cartesian(T& x, T& y, T& z, T r, T theta, T phi, T a)
{
  //
  // find the cartesian co-ordinates (aff,x,y,z) in the kerr spacetime
  //
  // result stored in first argument (T* x)
  //
  // arguments
  //   x      T*        stores the calcualted 'cartesian' co-ordinates
  //   r      T[4]      input Boyer-Lindquist co-ordinates
  //   a      T         BH spin parameter
  //
  x = sqrt(r*r + a*a)*sin(theta)*cos(phi);
  y = sqrt(r*r + a*a)*sin(theta)*sin(phi);
  z = r*cos(theta);
}

template<typename T>
T dot_product(T (*g)[4], T* u, T* v)
{
	//
	// returns the dot product of 4-vectors u and v under the metric g
	//
	T dotprod = 0;
	for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
		{
			dotprod += g[i][j] * u[i] * v[j];
		}

	return dotprod;
}

template<typename T>
void minkowski(T (*g)[4])
{
  //
  // return the Minkowski metric, eta = diag(1, -1, -1, -1)
  //
  // arguments:
  //   g      T* [4][4]      returned metric coefficients
  //
  for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
      g[i][j] = 0;

  g[0][0] = 1;
  g[1][1] = -1;
  g[2][2] = -1;
  g[3][3] = -1;
}

template<typename T>
void kerr_metric(T (*g)[4], T *r, T a)
{
  //
  // calculate Kerr metric coefficients at position r
  // (x0,x1,x2,x3) = (t,r,theta,phi)
  //
  // arguments:
  //   g      T* [4][4]      calculated metric coefficients
  //   r      T[4]           position
  //   a      T              BH spin parameter
  //
  const T rhosq = r[1]*r[1] + (a*cos(r[2]))*(a*cos(r[2]));
  const T delta = r[1]*r[1] - 2*r[1] + a*a;
  const T sigmasq = (r[1]*r[1] + a*a)*(r[1]*r[1] + a*a) - a*a*delta*sin(r[2])*sin(r[2]);

  // metric coefficients
  const T e2nu = rhosq * delta / sigmasq;
  const T e2psi = sigmasq * sin(r[2])*sin(r[2]) / rhosq;
  const T omega = 2*a*r[1] / sigmasq;

  for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
      g[i][j] = 0;

  g[0][0] = e2nu - omega*omega*e2psi;
  g[0][3] = omega*e2psi;
  g[3][0] = g[0][3];
  g[3][3] = -e2psi;
  g[1][1] = -rhosq/delta;
  g[2][2] = -rhosq;
}

template<typename T>
void tetrad(T* et, T* e1, T* e2, T* e3, T *r, T V, T a)
{
  //
  // calculate the tetrad basis vectors of an observer in the Kerr spacetime
  // at position r rotating at angular velocity (d(phi)/dt) V
  //
  // arguments:
  //   et      T* [4]      calculated time-like basis vector
  //   e1      T* [4]      calculated phi basis vector
  //   e2      T* [4]      calculated theta basis vector
  //   e3      T* [4]      calculated r basis vector
  //   r       T[4]        position of observer
  //   V       T           angular velocity (d(phi)/dt) of observer
  //   a       T           BH spin parameter
  //
  const T rhosq = r[1]*r[1] + (a*cos(r[2]))*(a*cos(r[2]));
  const T delta = r[1]*r[1] - 2*r[1] + a*a;
  const T sigmasq = (r[1]*r[1] + a*a)*(r[1]*r[1] + a*a) - a*a*delta*sin(r[2])*sin(r[2]);

  // metric coefficients
  const T e2nu = rhosq * delta / sigmasq;
  const T e2psi = sigmasq * sin(r[2])*sin(r[2]) / rhosq;
  const T omega = 2*a*r[1] / sigmasq;

  et[0] = (1/sqrt(e2nu))/sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu);
  et[1] = 0;
  et[2] = 0;
  et[3] = (1/sqrt(e2nu))*V / sqrt(1 - (V - omega)*(V - omega)*e2psi/e2nu);

  e1[0] = (V - omega)*sqrt(e2psi/e2nu) / sqrt(e2nu - (V-omega)*(V-omega)*e2psi);
  e1[1] = 0;
  e1[2] = 0;
  e1[3] = (1/sqrt(e2nu*e2psi))*(e2nu + V*omega*e2psi - omega*omega*e2psi) / sqrt(e2nu - (V-omega)*(V-omega)*e2psi);

  e2[0] = 0;
  e2[1] = 0;
  e2[2] = 1/sqrt(rhosq);
  e2[3] = 0;

  e3[0] = 0;
  e3[1] = sqrt(delta/rhosq);
  e3[2] = 0;
  e3[3] = 0;
}

template<typename T>
T lorentz(T* vel, T *v, T *r, T a)
{
  //
  // returns the Lorentz factor corresponding to a 4-velocity v
  // for a stationary observer at r
  //
  // arguments
  //   vel    T* [4]   stores 4-velocity as measured in observer's frame
  //   v      T[4]     input 4-velocity
  //   r      T[4]     position 4-vector of observer
  //   a      T        BH spin parameter
  //
  T g[4][4];
  T et[4], e1[4], e2[4], e3[4];
  T gvel[4];

  const T rhosq = r[1]*r[1] + (a*cos(r[2]))*(a*cos(r[2]));
  const T delta = r[1]*r[1] - 2*r[1] + a*a;
  const T sigmasq = (r[1]*r[1] + a*a)*(r[1]*r[1] + a*a) - a*a*delta*sin(r[2])*sin(r[2]);

  // observer rotating with frame-dragging
  const T omega = 2*a*r[1] / sigmasq;

  // metric coefficients
    kerr_metric(g, r, a);

  // tetrad basis vectors
    tetrad(et, e1, e2, e3, r, omega, a);

  // calculate components of v by dot product with basis vectors (gamma*c, gamma*vel)
  gvel[0] = g[0][0]*v[0]*et[0] + g[0][3]*v[0]*et[3] + g[3][0]*v[3]*et[0] + g[3][3]*v[3]*et[3];
  gvel[1] = g[0][0]*v[0]*e1[0] + g[0][3]*v[0]*e1[3] + g[3][0]*v[3]*e1[0] + g[3][3]*v[3]*e1[3];
  gvel[2] = g[2][2]*v[2]*e2[2];
  gvel[3] = g[1][1]*v[1]*e3[1];

  // calculate v from gamma*v
  for(int i=1; i<4; i++)
    vel[i] = gvel[i] / gvel[0];

  return gvel[0];   // v[0] is gamma (c=1)
}

template<typename T>
T disc_velocity_vector(T* v, T r, T a, int sign)
{
  //
  // calculate the 4-velocity of an element of the disc in stable circular orbit
  // in the equatorial plane of the Kerr geometry
  //
  // returns the angular velcocity (d(phi)/dt) measured in co-ordinate frame (T)
  //
  // arguments:
  //   v      T* [4]   calculated 4-velocity
  //   r      T        radius
  //   a      T        BH spin parameter
  //   sign   int           +1 for prograde or -1 for retrograde orbit
  //
  const T u = 1 / r;

  // constants of motion for stable circular orbit in Kerr spacetime
  const T k = ( 1 - 2*u + sign*a*sqrt(u*u*u) ) / sqrt( 1 - 3*u + sign*2*a*sqrt(u*u*u) );
  const T h = sign * ( 1 + a*a*u*u - sign*2*a*sqrt(u*u*u) ) / ( sqrt(u) * sqrt( 1 - 3*u + sign*2*a*sqrt(u*u*u) ) );

  v[0] = ( r*r*(r*r + a*a) + 2*a*a*r )*k - 2*a*r*h;
  v[0] = v[0] / ( r*r*(1 - (2/r))*(r*r + a*a) + 2*a*a*r);

  v[1] = 0;
  v[2] = 0;

  v[3] = 2*a*r*k + (r*r - 2*r)*h;
  v[3] = v[3] / ( r*r*(1-(2/r))*(r*r + a*a) + 2*a*a*r);

  // calculate and return angular velocity
  return v[3] / v[0];
}

template<typename T>
T disc_area(T r, T dr, T a)
{
  //
  // returns the proper area of an annulus in the equatorial plane
  // (T)
  //
  // arguments:
  //   r      T      radial co-ordinate of annulus
  //   dr     T      co-ordinate thickness of annulus
  //   a      T      BH spin parameter
  //
  const T rhosq = r*r;
  const T delta = r*r - 2*r + a*a;

  return sqrt( r*r + a*a + (2*a*a*r)/rhosq ) * sqrt(rhosq/delta)*dr;
}

template<typename T>
T rel_disc_area(T r, T dr, T a)
{
  //
  // return the effective area of a thin annulus in the equatorial plane (T)
  // taking into account relativistic effects on area
  //
  // arguments:
  //   r      T      r co-ordinate of annulus
  //   dr     T      thickness of annulus (<<r)
  //   a      T      BH spin parameter
  //
  T area;

  T g[4][4], et[4], e1[4], e2[4], e3[4], v[4], vel[3];
  T pos[] = {0,0,0,0};
  T V, gr_area, gamma;

  pos[1] = r;
    kerr_metric<T>(g, pos, a);
  V = disc_velocity_vector<T>(v, r, a, +1);
    tetrad(et, e1, e2, e3, pos, V, a);

  gr_area = disc_area<T>(r, dr, a);
  // Lorentz factor
  gamma = lorentz<T>(vel, v, pos, a);

  area = gr_area/gamma;

  return area;
}

template <typename T>
inline void momentum_from_consts(T& pt, T& pr, T& ptheta, T& pphi,
						T k, T h, T Q, int rdot_sign, int thetadot_sign,
						T r, T theta, T phi,
						const T a)
{
	//
	// Calculate photon momentum from constants of motion at a location
	//
	const T rhosq = r*r + (a*cos(theta))*(a*cos(theta));
	const T delta = r*r - 2*r + a*a;

	// tdot
	pt = (rhosq*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta))*k - 2*a*r*h;
	pt = pt / ( r*r * (1 + (a*cos(theta)/r)*(a*cos(theta)/r) - 2/r)*(r*r + a*a) + 2*a*a*r*sin(theta)*sin(theta) );

	// phidot
	pphi = 2*a*r*sin(theta)*sin(theta)*k + (r*r + (a*cos(theta))*(a*cos(theta)) - 2*r)*h;
	pphi = pphi / ( (r*r + a*a)*(r*r + (a*cos(theta))*(a*cos(theta)) - 2*r)*sin(theta)*sin(theta) + 2*a*a*r*sin(theta)*sin(theta)*sin(theta)*sin(theta) );

	// thetadot
	T thetadotsq = Q + (k*a*cos(theta) + h/tan(theta))*(k*a*cos(theta) - h/tan(theta));
	thetadotsq = thetadotsq / (rhosq*rhosq);

	// take the square roots and get the right signs
	ptheta = sqrt(abs(thetadotsq)) * thetadot_sign;

	// rdot
	T rdotsq = k*pt - h*pphi - rhosq*ptheta*ptheta;
	rdotsq = rdotsq * delta/rhosq;

	pr = sqrt(abs(rdotsq)) * rdot_sign;
}


#endif /* KERR_H_ */
