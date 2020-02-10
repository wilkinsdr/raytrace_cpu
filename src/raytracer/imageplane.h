/*
 * imagePlane.h
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#ifndef IMAGEPLANE_H_
#define IMAGEPLANE_H_

#include "raytracer.h"

template <typename T>
class ImagePlane : public Raytracer<T>
{
private:
	T D;
	T incl;
	T phi0;

	int Nx, Ny;

	T* m_plane_x;
	T* m_plane_y;

public:
	ImagePlane( T dist, T inc, T x0, T xmax, T dx, T y0, T ymax, T dy, T spin, T tol, T phi = 0);
	~ImagePlane();

	void init_image_plane(T D, T incl, T phi0, T x0, T xmax, T dx, T y0, T ymax, T dy);

	void redshift_start( );
	void redshift(bool projradius);

};

#endif /* IMAGEPLANE_H_ */
