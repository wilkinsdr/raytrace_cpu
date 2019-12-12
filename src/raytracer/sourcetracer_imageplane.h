/*
 * imagePlane.h
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#ifndef SOURCETRACER_IMAGEPLANE_H_
#define SOURCETRACER_IMAGEPLANE_H_

#include "source_tracer.h"

template <typename T>
class SourceTracer_ImagePlane : public SourceTracer<T>
{
private:
	T D;
	T incl;
	T phi0;

	int Nx, Ny;

	T* m_plane_x;
	T* m_plane_y;

public:
	SourceTracer_ImagePlane( T dist, T inc, T x0, T xmax, T dx, T y0, T ymax, T dy, T spin, T init_en0, T init_enmax, int init_Nen, bool init_logbin_en = false, T tol = TOL, T phi = 0 );
	~SourceTracer_ImagePlane();

	void InitImagePlane( T D, T incl, T phi0, T x0, T xmax, T dx, T y0, T ymax, T dy);

	void RedshiftStart( );
	void Redshift( bool projradius );

};

#endif /* SOURCETRACER_IMAGEPLANE_H_ */
