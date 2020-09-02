/*
 * imagePlane.h
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#ifndef Mapper_IMAGEPLANE_H_
#define Mapper_IMAGEPLANE_H_

#include "mapper.h"

class Mapper_ImagePlane : public Mapper
{
private:
	double D;
	double incl;
	double phi0;

	int Nx, Ny;

	double* m_plane_x;
	double* m_plane_y;

	virtual long get_num_rays()
    {
	    return Nx*Ny;
    }

public:
    Mapper_ImagePlane( double dist, double inc, double x0, double xmax, double dx, double y0, double ymax, double dy, double spin, double init_r0, double init_rmax, int init_Nr, int init_Ntheta, int init_Nphi, bool init_logbin_r, double init_thetamax = M_PI_2, double tol = TOL, double phi = 0 );
	~Mapper_ImagePlane();

	void init_image_plane( double D, double incl, double phi0, double x0, double xmax, double dx, double y0, double ymax, double dy);

	void redshift_start( );
	void redshift( bool projradius );

};

#endif /* Mapper_IMAGEPLANE_H_ */
