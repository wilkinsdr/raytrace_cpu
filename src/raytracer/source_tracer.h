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

#ifndef RAYTRACER_H_
#define RAYTRACER_H_

// set a default tolerance (precision used for variable integration step sizes)
#define TOL 100
// set a default outer radius for ray propagation
#define RLIM 1000

// total number of steps allowed per ray before it's aborted
#define STEPLIM 10000000
// smallest integration step  to prevent infinite loops of infinitesimal steps
#define MIN_STEP 1E-12
// minimum number of integration steps needed for a ray to be included in analysis
#define COUNT_MIN 100

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include "../include/kerr.h"
#include "../include/text_output.h"
#include "raytracer.h"


template <typename T>
class SourceTracer : public RayTracer
{
private:
	float tolerance;
	float horizon;

protected:	// these members need to be accessible by derived classes to set up different X-ray sources
	T **emis, **absorb;

public:
    SourceTracer( int num_rays, float spin, float tol = TOL);
    ~SourceTracer( );

    void run_source_trace( T r_max = 1000, T theta_max = M_PI/2, TextOutput* outfile = 0, int write_step = 1, T write_rmax = -1, T write_rmin = -1, bool write_cartesian = true );
    inline int propagate_source(int ray, const T rlim, const T thetalim, const int steplim, TextOutput* outfile = 0, int write_step = 1, T write_rmax = -1, T write_rmin = -1, bool write_cartesian = true);

};

#endif /* RAYTRACER_H_ */
