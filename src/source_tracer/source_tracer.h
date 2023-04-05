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
#include <functional>
using namespace std;

#include "../include/kerr.h"
#include "../include/text_output.h"
#include "../include/progress_bar.h"
#include "../raytracer/raytracer.h"
//#include "mapper.h"


template <typename T>
class SourceTracer : public Raytracer<T>
{
private:
	T source_size_xy, source_size_z, source_vel;
	int source_motion;

	bool reverse;

	bool stopping_fn_set;
	//bool (*stopping_fn)(T, T, T, T, T*);
	//function<bool(T, T, T, T, T*)> stopping_fn;
	bool (*stopping_fn)(T, T, T, T, T*);
	T *stopping_args;

	//Mapper *mapper;

public:
    SourceTracer( int num_rays, float spin_par, T init_en0, T init_enmax, int init_Nen, bool init_logbin_en = false, T init_t0 = 0, T init_tmax = 1E6, int init_Nt = 1, float toler = TOL, bool reverse = false );
    ~SourceTracer( );

	T *energy;
	T en0, enmax, den, t0, tmax, dt;
	int Nen, Nt;
	bool logbin_en;

	T **emis, **absorb, **emis_ent;

    void run_source_trace( T r_max = 1000, T theta_max = M_PI_2, int show_progress = 1, TextOutput* outfile = 0, int write_step = 1, T write_rmax = -1, T write_rmin = -1, bool write_cartesian = true );
    inline int propagate_source(int ray, const T rlim, const T thetalim, const int steplim, TextOutput* outfile = 0, int write_step = 1, T write_rmax = -1, T write_rmin = -1, bool write_cartesian = true);

	void set_source(T size_xy, T size_z, T vel, int motion = 1)
	{
		source_size_xy = size_xy;
		source_size_z = size_z;
		source_motion = motion;
		source_vel = vel;
	}

//	void set_stopping_fn(function<bool(T, T, T, T, T*)> fn, T* args)
//	{
//		stopping_fn = fn;
//		stopping_fn_set = true;
//		stopping_args = args;
//	}

    void set_stopping_fn(bool (*fn)(T, T, T, T, T*), T* args)
    {
        stopping_fn = fn;
        stopping_fn_set = true;
        stopping_args = args;
    }

//	void add_mapper(Mapper* map_ptr)
//	{
//		mapper = map_ptr;
//	}
};

#endif /* SOURCETRACER_H_ */
