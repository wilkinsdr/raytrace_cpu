/*
 * imagePlane.h
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#ifndef Mapper_IMAGEPLANE_H_
#define Mapper_IMAGEPLANE_H_

#include "mapper.h"

template <typename T>
class Mapper_ImagePlane : public Mapper<T>
{
private:
	T D;
	T incl;
	T phi0;

	int Nx, Ny;

	T* m_plane_x;
	T* m_plane_y;

public:
	Mapper_ImagePlane( T dist, T inc, T x0, T xmax, T dx, T y0, T ymax, T dy, T spin, T init_en0, T init_enmax, int init_Nen, bool init_logbin_en = false, T tol = TOL, T phi = 0 );
	~Mapper_ImagePlane();

	void init_image_plane( T D, T incl, T phi0, T x0, T xmax, T dx, T y0, T ymax, T dy);

	void redshift_start( );
	void redshift( bool projradius );

};

#endif /* Mapper_IMAGEPLANE_H_ */
