/*
 * pointsource_vel.h
 *
 *  Created on: 13 Jul 2017
 *      Author: drw
 */

#ifndef POINTSOURCE_VEL_H_
#define POINTSOURCE_VEL_H_

#include "raytracer.h"

template <typename T>
class PointSourceVel : public Raytracer<T>
{
private:
	T energy;
	T velocity;

	int n_cosalpha;
	int n_beta;

	T* m_cosalpha;
	T* m_beta;

    T m_tetrad[4][4];

public:
	PointSourceVel( T* pos, T* V, T spin, T tol, T dcosalpha, T dbeta, T cosalpha0 = -0.999999, T cosalphamax = 0.995, T beta0 = -0.995*M_PI, T betamax = M_PI, T E = 1 );
	~PointSourceVel();

	void InitPointSourceVel( T* pos, T dcosalpha, T dbeta, T cosalpha0 = -0.999999, T cosalphamax = 0.995, T beta0 = -0.995*M_PI, T betamax = M_PI );

	void SolveTetrad( T* pos, T* V, T spin );

	void RedshiftStart( );
	void Redshift( T V );

};

#endif /* POINTSOURCE_VEL_H_ */
