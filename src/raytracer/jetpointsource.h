/*
 * pointsource.h
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#ifndef JETPOINTSOURCE_H_
#define JETPOINTSOURCE_H_

#include "raytracer.h"

template <typename T>
class JetPointSource : public Raytracer<T>
{
private:
	T energy;
	T velocity;

	int n_cosalpha;
	int n_beta;

	T* m_cosalpha;
	T* m_beta;

	inline void CalculateConstants_RadialMotion(int ray, T alpha, T beta, T V, T E);

public:
	JetPointSource( T* pos, T V, T spin, T tol, T dcosalpha, T dbeta, T cosalpha0 = -0.999999, T cosalphamax = 0.995, T beta0 = -0.995*M_PI, T betamax = M_PI, T E = 1 );
	~JetPointSource( );

	void init_jet_pointsource(T* pos, T dcosalpha, T dbeta, T cosalpha0 = -0.999999, T cosalphamax = 0.995, T beta0 =
    -0.995 * M_PI, T betamax = M_PI);

	void redshift_start( );
	void redshift(T V);

};

#endif /* JETPOINTSOURCE_H_ */
