//
// Created by Will Surgent on 8/8/23.
//

#ifndef RAYTRACE_CPU_LOG_IMAGEPLANE_H
#define RAYTRACE_CPU_LOG_IMAGEPLANE_H

#include "raytracer.h"

template <typename T>
class LogImagePlane : public Raytracer<T>
{
private:
    T D;
    T incl;
    T phi0;

    T m_x0, m_xmax, m_dx;
    T m_y0, m_ymax, m_dy;

    int Nx, Ny;
    int m_nx, m_ny;
    int m_quad;

public:
    //LogImagePlane(T dist, T inc, T x0, T xmax, T dx, T y0, T ymax, T dy, T spin, T quad, T phi0, T excise, T precision = PRECISION, T max_tstep);
    LogImagePlane(T dist, T inc, T x0, T xmax, T dx, T y0, T ymax, T dy, T spin, int quad = 0, T phi0 = 0, T excise = 0, T precision = PRECISION, T max_tstep = MAXDT);
//	~ImagePlane();

//    void ray_sort(const int threads, T rays[]);

//    void init_image_plane(T D, T incl, T phi0, T x0, T xmax, T dx, T y0, T ymax, T dy);

    void init_log_imageplane(T dist, T incl, T x0, T xmax, T dx, T y0, T ymax, T dy, T spin, int quad = 0, T phi0 = 0,
                       T excise = 0, T precision = PRECISION);

    void redshift_start( );
    void redshift(RayDestination<T>* destination, bool projradius);

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
#endif //RAYTRACE_CPU_LOG_IMAGEPLANE_H
