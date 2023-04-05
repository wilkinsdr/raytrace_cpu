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

    T m_x0, m_xmax, m_dx;
    T m_y0, m_ymax, m_dy;

	int Nx, Ny;

public:
	ImagePlane( T dist, T inc, T x0, T xmax, T dx, T y0, T ymax, T dy, T spin, T phi, T precision = PRECISION);
//	~ImagePlane();

	void init_image_plane(T D, T incl, T phi0, T x0, T xmax, T dx, T y0, T ymax, T dy);

	void redshift_start( );
	void redshift(bool projradius);

    inline int get_x_index(int ix)
    {
        //
        // returns the index along the X direction of the 2D image place given the index of the 1D array element in raytracer variables
        //
        return static_cast<int>( ix / Ny );
    }

    inline int get_y_index(int ix)
    {
        //
        // returns the index along the Y direction of the 2D image place given the index of the 1D array element in raytracer variables
        //
        return ix % Ny;
    }

    inline T ray_x(int ix)
    {
        return m_x0 + get_x_index(ix) * m_dx;
    }

    inline T ray_y(int ix)
    {
        return m_y0 + get_y_index(ix) * m_dy;
    }
};

#endif /* IMAGEPLANE_H_ */
