/*
 * pointsource.h
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#ifndef HEALPIX_POINTSOURCE_H_
#define HEALPIX_POINTSOURCE_H_

#include "raytracer.h"

#include "../include/healpix.h"

template <typename T>
class HealpixPointSource : public Raytracer<T>
{
private:
	T energy;
	T velocity;

    int order;
    int nside, npix;

public:
	HealpixPointSource( T* pos, T V, T spin, int order, int motion = 0, int basis = 0, T tol = TOL );

	void init_healpix_pointsource(T* pos, int order, int motion = 0, int basis = 0 );

	void redshift_start( );
	void redshift(T V );

	int get_num_pix()
    {
	    return npix;
    }

};

#endif /* HEALPIX_POINTSOURCE_H_ */
