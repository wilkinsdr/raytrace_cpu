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

#ifndef SOURCETRACER_H_
#define SOURCETRACER_H_

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include "../include/kerr.h"
#include "../include/text_output.h"
#include "raytracer.h"


template <typename T>
class SourceTracer : public Raytracer<T>
{
private:
	T source_size_xy, source_size_z, source_vel;
	int source_motion;

	bool reverse;

public:
    SourceTracer( int num_rays, float spin_par, T init_en0, T init_enmax, int init_Nen, bool init_logbin_en = false, float toler = TOL, bool reverse = false );
    ~SourceTracer( );

	T *energy;
	T en0, enmax, den;
	int Nen;
	bool logbin_en;

	T **emis, **absorb;

    void run_source_trace( T r_max = 1000, T theta_max = M_PI/2, TextOutput* outfile = 0, int write_step = 1, T write_rmax = -1, T write_rmin = -1, bool write_cartesian = true );
    inline int propagate_source(int ray, const T rlim, const T thetalim, const int steplim, TextOutput* outfile = 0, int write_step = 1, T write_rmax = -1, T write_rmin = -1, bool write_cartesian = true);

	void set_source(T size_xy, T size_z, T vel, int motion = 1)
	{
		source_size_xy = size_xy;
		source_size_z = size_z;
		source_motion = motion;
		source_vel = vel;
	}
};

#endif /* SOURCETRACER_H_ */
