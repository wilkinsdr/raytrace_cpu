/*
 * pointsource.h
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#ifndef POINTSOURCE_H_
#define POINTSOURCE_H_

#include "raytracer.h"

template <typename T>
class Beam : public Raytracer<T>
{
private:
	T energy;
	T velocity;

	int n_cosalpha;
	int n_beta;

	T* m_cosalpha;
	T* m_beta;

public:
	PointSource( T* pos, T V, T spin, T tol, T dcosalpha, T dbeta, T cosalpha0 = -0.999999, T cosalphamax = 0.995, T beta0 = -0.995*M_PI, T betamax = M_PI, T E = 1 );
	~PointSource();

	void InitPointSource( T* pos, T dcosalpha, T dbeta, T cosalpha0 = -0.999999, T cosalphamax = 0.995, T beta0 = -0.995*M_PI, T betamax = M_PI );

	void RedshiftStart( );
	void Redshift( T V );

};

#endif /* POINTSOURCE_H_ */
