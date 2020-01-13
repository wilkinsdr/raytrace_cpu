/*
 * pointsource.h
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#ifndef MAPPER_POINTSOURCE_H_
#define MAPPER_POINTSOURCE_H_

#include "mapper.h"

template <typename T>
class Mapper_PointSource : public Mapper<T>
{
private:
	T energy;
	T velocity;

	int n_cosalpha;
	int n_beta;

	T* m_cosalpha;
	T* m_beta;

public:
	Mapper_PointSource( T* pos, T V, T spin, T dcosalpha, T dbeta, T init_r0, T init_rmax, int init_Nr, int init_Ntheta, int init_Nphi, bool init_logbin_r, T init_thetamax, T cosalpha0 = -0.999999, T cosalphamax = 0.995, T beta0 = -0.995*M_PI, T betamax = M_PI, T tol = TOL, T E = 1 );
	~Mapper_PointSource();

	void InitPointSource( T* pos, T dcosalpha, T dbeta, T cosalpha0 = -0.999999, T cosalphamax = 0.995, T beta0 = -0.995*M_PI, T betamax = M_PI );

	void RedshiftStart( );
	void Redshift( T V );

};

#endif /* MAPPER_POINTSOURCE_H_ */
