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
//#define STEPLIM 1000
// smallest integration step  to prevent infinite loops of infinitesimal steps
#define MIN_STEP 1E-6
// minimum number of integration steps needed for a ray to be included in analysis
#define COUNT_MIN 10000

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include "../include/kerr.h"
#include "../include/text_output.h"
#include "../include/progress_bar.h"

template <typename T>
class Raytracer
{
protected:	// these members need to be accessible by derived classes to set up different X-ray sources
	int nRays;
	T spin;
	float tolerance;
	float horizon;
	float max_phistep;
	float max_tstep;

	// pointers for ray variables
	T *m_t, *m_r, *m_theta, *m_phi;
	T *m_pt, *m_pr, *m_ptheta, *m_pphi;
	T *m_k, *m_h, *m_Q;
	T *m_emit, *m_redshift;
	int *m_steps;
	int *m_status;
	int *m_rdot_sign, *m_thetadot_sign;

	inline void CalculateConstants(int ray, T alpha, T beta, T V, T E);
	inline void CalculateConstantsFromP(int ray, T pt, T pr, T ptheta, T pphi);

public:
    Raytracer( int num_rays, float spin, float tol = TOL, float init_max_phistep = 0.1, float init_max_tstep = 1 );
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

    T* get_time( )
    {
    	return m_t;
    }

    T* get_r( )
    {
    	return m_r;
    }

    T* get_theta( )
    {
    	return m_theta;
    }

    T* get_phi( )
    {
    	return m_phi;
    }

    int* get_steps( )
    {
    	return m_steps;
    }

    T* get_redshift( )
    {
    	return m_redshift;
    }
	int* get_status( )
	{
		return m_status;
	}

    void map_results(int*& p_steps, T*& p_t, T*& p_r, T*& p_theta, T*& p_phi, T*& p_redshift )
    {
    	//
    	// Maps the pointers to the raytrace variable arrays for access by external code
    	//
    	p_steps = m_steps;
    	p_t = m_t;
    	p_r = m_r;
    	p_theta = m_theta;
    	p_phi = m_phi;
    	p_redshift = m_redshift;
    }
    void map_momentum(T*& p_pt, T*& p_pr, T*& p_ptheta, T*& p_pphi )
    {
    	//
    	// Maps the pointers to the momentum variable arrays for access by external code
    	//
    	p_pt = m_pt;
    	p_pr = m_pr;
    	p_ptheta = m_ptheta;
    	p_pphi = m_pphi;
    }
    void map_consts(T*& p_k, T*& p_Q, T*& p_h )
    {
    	//
    	// Maps the pointers to the arrays containing the constants of motion for access by external code
    	//
    	p_k = m_k;
    	p_Q = m_Q;
    	p_h = m_h;
    }

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

};

#endif /* RAYTRACER_H_ */
