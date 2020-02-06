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
	HealpixPointSource( T* pos, T V, T spin, int order, T tol = TOL );

	void InitHealpixPointSource( T* pos, int order );

	void RedshiftStart( );
	void Redshift( T V );

	int GetNumPix()
    {
	    return npix;
    }

};

#endif /* HEALPIX_POINTSOURCE_H_ */
