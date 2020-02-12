//
// Created by drw on 01/03/18.
//

#ifndef CUDAKERR_DISC_H
#define CUDAKERR_DISC_H

#include "kerr.h"
#include "gramschmidt_basis.h"

template <typename T>
T rel_vector_disc_area(T r, T dr, T dphi, T a)
{
	T et[4], e1[4], e2[4], e3[4];
	T g[4][4];
	T pos[] = {0, r, M_PI/2, 0};

	T V = DiscVelocity(r, a, +1);

	Tetrad(et, e1, e2, e3, pos, V, a);
	Metric(g, pos, a);

	T side1[] = {0, dr, 0, 0};
	T side2[] = {0, 0, 0, dphi};

	T side1_prime[] = { DotProduct(g, side1, et), DotProduct(g, side1, e1), DotProduct(g, side1, e2), DotProduct(g, side1, e3)};
	T side2_prime[] = { DotProduct(g, side2, et), DotProduct(g, side2, e1), DotProduct(g, side2, e2), DotProduct(g, side2, e3)};

	T cross[] = {side1_prime[2]*side2_prime[3] - side1_prime[3]*side2_prime[2], side1_prime[3]*side2_prime[1] - side1_prime[1]*side2_prime[3], side1_prime[1]*side2_prime[2] - side1_prime[2]*side2_prime[1] };
	// return the full parallelogram area
	return sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
}

template <typename T>
T rel_vector_disc_area_plunge(T r, T dr, T dphi, T a)
{
	T *et, *e1, *e2, *e3;
	T g[4][4];
	T pos[] = {0, r, M_PI/2, 0};
	T vel[4];

	const T delta = r*r - 2*r + a*a;

	const T r_isco = KerrISCO(a, +1);
	// constants of motion for the ISCO
	const T u = 1/r_isco;
	const T k = (1 - 2*u + a*u*sqrtf(u)) / sqrtf(1 - 3*u + 2*a*u*sqrtf(u));
	const T h = (1 + a*a*u*u - 2*a*u*sqrtf(u)) / sqrtf( u*(1 - 3*u + 2*a*u*sqrtf(u)) );

	// 4-velocity of a plunging orbit from the ISCO
	vel[0] = (1/delta) * ( (r*r + a*a + 2*a*a/r)*k - 2*a*h/r );
	vel[1] = -1*sqrtf( k*k - 1 + 2/r + (a*a*(k*k - 1) - h*h)/(r*r) + 2*(h-a*k)*(h-a*k)/(r*r*r) );
	vel[2] = 0;
	vel[3] = (1/delta) * (2*a*k/r + (1 - 2/r)*h);

	// if numerical error pushes v_r sqrt operand below zero, we're at the ISCO so v_r = 0
	if(!(abs(vel[1]) > 0)) vel[1] = 0;

	GramSchmidt_Basis<double> basis(pos, vel, a);
	et = basis.vectors[0];
	e1 = basis.vectors[1];
	e2 = basis.vectors[2];
	e3 = basis.vectors[3];

	Metric(g, pos, a);

	T side1[] = {0, dr, 0, 0};
	T side2[] = {0, 0, 0, dphi};

	T side1_prime[] = { DotProduct(g, side1, et), DotProduct(g, side1, e1), DotProduct(g, side1, e2), DotProduct(g, side1, e3)};
	T side2_prime[] = { DotProduct(g, side2, et), DotProduct(g, side2, e1), DotProduct(g, side2, e2), DotProduct(g, side2, e3)};

	T cross[] = {side1_prime[2]*side2_prime[3] - side1_prime[3]*side2_prime[2], side1_prime[3]*side2_prime[1] - side1_prime[1]*side2_prime[3], side1_prime[1]*side2_prime[2] - side1_prime[2]*side2_prime[1] };
	// return the full parallelogram area
	return sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
}

template <typename T>
T rel_vector_disc_area_plunge_varradius(T r, T dr, T dphi, T a, T r_plunge = -1)
{
	T *et, *e1, *e2, *e3;
	T g[4][4];
	T pos[] = {0, r, M_PI/2, 0};
	T vel[4];

	const T delta = r*r - 2*r + a*a;

	const T r_isco = KerrISCO(a, +1);

	if(r_plunge < 0) r_plunge = r_isco;

	// constants of motion for the ISCO
	const T u = 1/r_plunge;
	const T k = (1 - 2*u + a*u*sqrtf(u)) / sqrtf(1 - 3*u + 2*a*u*sqrtf(u));
	const T h = (1 + a*a*u*u - 2*a*u*sqrtf(u)) / sqrtf( u*(1 - 3*u + 2*a*u*sqrtf(u)) );

	// 4-velocity of a plunging orbit from the ISCO
	vel[0] = (1/delta) * ( (r*r + a*a + 2*a*a/r)*k - 2*a*h/r );
	vel[1] = -1*sqrtf( k*k - 1 + 2/r + (a*a*(k*k - 1) - h*h)/(r*r) + 2*(h-a*k)*(h-a*k)/(r*r*r) );
	vel[2] = 0;
	vel[3] = (1/delta) * (2*a*k/r + (1 - 2/r)*h);

	// if numerical error pushes v_r sqrt operand below zero, we're at the ISCO so v_r = 0
	if(!(abs(vel[1]) > 0)) vel[1] = 0;

	GramSchmidt_Basis<double> basis(pos, vel, a);
	et = basis.vectors[0];
	e1 = basis.vectors[1];
	e2 = basis.vectors[2];
	e3 = basis.vectors[3];

	Metric(g, pos, a);

	T side1[] = {0, dr, 0, 0};
	T side2[] = {0, 0, 0, dphi};

	T side1_prime[] = { DotProduct(g, side1, et), DotProduct(g, side1, e1), DotProduct(g, side1, e2), DotProduct(g, side1, e3)};
	T side2_prime[] = { DotProduct(g, side2, et), DotProduct(g, side2, e1), DotProduct(g, side2, e2), DotProduct(g, side2, e3)};

	T cross[] = {side1_prime[2]*side2_prime[3] - side1_prime[3]*side2_prime[2], side1_prime[3]*side2_prime[1] - side1_prime[1]*side2_prime[3], side1_prime[1]*side2_prime[2] - side1_prime[2]*side2_prime[1] };
	// return the full parallelogram area
	return sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
}

template <typename T>
T integrate_disc_area(T rmin, T rmax, T a, bool force_keplerian = false, int Nr = 50, T dphi = 0.1, bool logbin_r = true)
{
	const T dr = (logbin_r) ? exp(log(rmax/rmin)/(Nr-1)) : (rmax - rmin)/(Nr - 1);

	const T r_isco = KerrISCO(a, +1);

	T area = 0;
	for(T r=rmin; r<rmax; r = (logbin_r) ? (r*dr) : (r+dr))
	{
		const T this_dr = (logbin_r) ? r*(dr-1) : dr;
		const T this_area = (r >= r_isco || force_keplerian) ? rel_vector_disc_area(r, this_dr, dphi, a) : rel_vector_disc_area_plunge(r, this_dr, dphi, a);
		if(this_area > 0) area += this_area;
	}

	return area;
}

template <typename T>
T integrate_disc_area_varplungeradius(T rmin, T rmax, T a, T r_plunge = -1, int Nr = 50, T dphi = 0.1, bool logbin_r = true)
{
	const T dr = (logbin_r) ? exp(log(rmax/rmin)/(Nr-1)) : (rmax - rmin)/(Nr - 1);

	const T r_isco = KerrISCO(a, +1);

	if(r_plunge < 0) r_plunge = r_isco;

	T area = 0;
	for(T r=rmin; r<rmax; r = (logbin_r) ? (r*dr) : (r+dr))
	{
		const T this_dr = (logbin_r) ? r*(dr-1) : dr;
		const T this_area = (r >= r_plunge) ? rel_vector_disc_area(r, this_dr, dphi, a) : rel_vector_disc_area_plunge_varradius(r, this_dr, dphi, a, r_plunge);
		if(this_area > 0) area += this_area;
	}

	return area;
}

#endif //CUDAKERR_DISC_H
