/*
 * pointsource.h
 *
 *  Created on: 25 Apr 2016
 *      Author: drw
 */

#ifndef MAPPER_POINTSOURCE_H_
#define MAPPER_POINTSOURCE_H_

#include "mapper.h"

class Mapper_PointSource : public Mapper
{
private:
	double energy;
	double velocity;

	int n_cosalpha;
	int n_beta;

	double* m_cosalpha;
	double* m_beta;

public:
	Mapper_PointSource( double* pos, double V, double spin, double dcosalpha, double dbeta, double init_r0, double init_rmax, int init_Nr, int init_Ntheta, int init_Nphi, bool init_logbin_r, double init_thetamax, double cosalpha0 = -0.999999, double cosalphamax = 0.995, double beta0 = -0.995*M_PI, double betamax = M_PI, double tol = TOL, double E = 1 );
	~Mapper_PointSource();

	void InitPointSource( double* pos, double dcosalpha, double dbeta, double cosalpha0 = -0.999999, double cosalphamax = 0.995, double beta0 = -0.995*M_PI, double betamax = M_PI );

	void RedshiftStart( );
	void Redshift( double V );

};

#endif /* MAPPER_POINTSOURCE_H_ */
