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

// set a default precision used for variable integration step sizes based on distance from event horizon
#define PRECISION 100
#define TOL 100
//
#define THETA_PRECISION 50
// default maximum allowed step in co-ordinate time
#define MAXDT 1
// only obey this inside this radius so rays don't take ages a long way from the black hole
#define MAXDT_RLIM 100
// default maximum allowed step in phi
#define MAXDPHI 0.1
// set a default outer radius for ray propagation
#define RLIM 1000
// total number of steps allowed per ray before it's aborted
#define STEPLIM 10000000
// number of integration steps per GPU thread before integration is paused and kernel must be called again
// to avoid thread time limits on GPUs running X servers
#define THREAD_STEPLIM 10000000
// smallest integration step  to prevent infinite loops of infinitesimal steps
#define MIN_STEP 1E-3
// minimum number of integration steps before sign of rdot and thetadot are allowed to change
#define COUNT_MIN 100

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include "../include/kerr.h"
#include "../include/text_output.h"
#include "../include/progress_bar.h"

template <typename T>
struct Ray
{
    T t, r, theta, phi;
    T pt, pr, ptheta, pphi;
    T k, h, Q;
    T emit, redshift;
    int steps;
    int status;
    int rdot_sign, thetadot_sign;
    T alpha, beta;
};

template <typename T>
class Raytracer
{
private:
    T precision;
    T theta_precision;
    T max_tstep;
    T max_phistep;
    T maxtstep_rlim;

protected:	// these members need to be accessible by derived classes to set up different X-ray sources
	int nRays;
	T spin;
	T horizon;

	inline void calculate_constants(int ray, T alpha, T beta, T V, T E);
	inline void calculate_constants_from_p(int ray, T pt, T pr, T ptheta, T pphi);

public:
    Ray<T> *rays;

    Raytracer( int num_rays, T spin, T precision = PRECISION, T init_max_phistep = MAXDPHI, T init_max_tstep = MAXDT );
    ~Raytracer( );

    void run_raytrace(T r_max = 1000, T theta_max = M_PI / 2, int show_progress = 1, TextOutput* outfile = 0
                      , int write_step = 1, T write_rmax = -1, T write_rmin = -1, bool write_cartesian = true);
    inline int propagate(int ray, const T rlim, const T thetalim, const int steplim, TextOutput* outfile = 0
                         , int write_step = 1, T write_rmax = -1, T write_rmin = -1, bool write_cartesian = true);

    void redshift_start(T V, bool reverse = false, bool projradius = false);
    void redshift(T V, bool reverse = false, bool projradius = false, int motion = 0);
	inline T ray_redshift( T V, bool reverse, bool projradius, T r, T theta, T phi, T k, T h, T Q, int rdot_sign, int thetadot_sign, T emit, int motion = 0 );

    void range_phi(T min = -1 * M_PI, T max = M_PI);

    void calculate_momentum( );

    int get_count( )
    {
    	//
    	// Returns the total number of threads (number of entries in raytrace variable arrays).
    	// This is not necessarily the number of rays as those out of the specified bounds will not have been traced,
    	// though this should be used for iterating over the raytrace results.
    	//
    	return nRays;
    }

    void set_boundary(T r = -1)
    {
    	//
		// Sets a hard boundary (inner radius) through which rays cannot pass, e.g. for raytracing around neutron stars
    	// If called with no argument, sets the boundary to the event horizon.
		//
    	if(r > 0)
    		horizon = r;
    	else
    		horizon = kerr_horizon<T>(spin);
    }

    T calculate_horizon( )
    {
    	return kerr_horizon<T>(spin);
    }

    void set_precision(T precision)
    {
        precision = precision;
    }

    void set_precision(T precision = PRECISION, T theta_precision = THETA_PRECISION)
    {
        precision = precision;
        theta_precision = theta_precision;
    }

    void set_max_tstep(T max, T rlim = MAXDT_RLIM)
    {
        max_tstep = max;
        maxtstep_rlim = rlim;
    }

    void set_max_phistep(T max)
    {
        max_phistep = max;
    }

};

#endif /* RAYTRACER_H_ */
