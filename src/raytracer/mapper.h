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
    Mapper(int num_rays, float spin_par, T init_r0, T init_rmax, int init_Nr, int init_Ntheta, int init_Nphi, bool init_logbin_r, T init_thetamax = M_PI_2, float toler = TOL, bool reverse = false);
    Mapper(char* load_filename);
    ~Mapper( );

	T r0, rmax, theta_max;
	int Nr, Ntheta, Nphi;
	T dr, dtheta, dphi;
	bool logbin_r;

	Array3D<T> *map_time, *map_redshift, *map_flux;
	Array3D<int> *map_Nrays;

    void run_map( T r_max = 1000 );
    inline int map_ray(int ray, const T rlim, const T thetalim, const int steplim);
    void average_rays();

    void save(char* filename);

};

#endif /* MAPPER_H_ */
