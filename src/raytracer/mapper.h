/*
 * gpu_raytracer.h
 *
 *  Created on: 4 Sep 2013
 *      Author: drw
 *
 *  General relativistic ray tracing in the Kerr spacetime
 *
 *  This class provides the basic functionality fundamental to the rays (e.g. propagation, redshift calculation)
 *  and is used as the base class for specific types of source for ray tracing in different scenarios.
 *
 *  This class is a template to set the floating point precision (e.g. float or double)
 */

#ifndef MAPPER_H_
#define MAPPER_H_

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include "../include/kerr.h"
#include "../include/text_output.h"
#include "../include/array.h"
#include "raytracer.h"


template <typename T>
class Mapper : public Raytracer<T>
{
private:
	int motion;
	int vel_mode;
	T vel;

	bool reverse;



public:
    Mapper( int num_rays, float spin_par, T init_en0, T init_enmax, int init_Nen, bool init_logbin_en = false, float toler = TOL, bool reverse = false );
    ~Mapper( );

	T r0, rmax;
	int Nr, Ntheta, Nphi;

	Array3D<T> *map_t, *map_redshift, *map_flux;
	Array3D<int> *map_Nrays;

    void run_map( T r_max = 1000, T theta_max = M_PI/2, TextOutput* outfile = 0, int write_step = 1, T write_rmax = -1, T write_rmin = -1, bool write_cartesian = true );
    inline int map_ray(int ray, const T rlim, const T thetalim, const int steplim, TextOutput* outfile = 0, int write_step = 1, T write_rmax = -1, T write_rmin = -1, bool write_cartesian = true);

};

#endif /* MAPPER_H_ */
