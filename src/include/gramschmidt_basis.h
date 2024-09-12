//
// Created by drw on 01/03/18.
//

#ifndef CUDAKERR_GRAMSCHMIDT_BASIS_H
#define CUDAKERR_GRAMSCHMIDT_BASIS_H

#include "kerr.h"
#include <cmath>

template <typename T>
class BasisVectors
{
public:
	T vectors[4][4];
};

template <typename T>
class GramSchmidt_Basis : public BasisVectors<T>
{
private:
	T g[4][4];

public:
	GramSchmidt_Basis(T* pos, T* vel, T spin)
	{
		kerr_metric(g, pos, spin);
		orthogonalise(BasisVectors<T>::vectors, pos, vel);
		normalise();
	}

	void orthogonalise(T (*vec)[4], T* pos, T* vel)
	{
		// v is the initial giess
		double v[4][4];
		// e is the working set of basis vectors as we orthogonalise them
		double e[4][4];

		// initial guess of the basis vectors
		// v0 is the 4-velocity
		v[0][0] = vel[0];
		v[0][1] = vel[1];
		v[0][2] = vel[2];
		v[0][3] = vel[3];
		// 1, 2 and 3 are initially in the phi, theta and r directions
		v[1][1] = 1; // put the radial one first to prioritise getting one in approximately this direction
		v[2][2] = 1;
		v[3][3] = 1;

		// u0 = v0
		for(int a=0; a<4; a++)
			e[0][a] = v[0][a];

		// Gram-Schmidt orthogonalisation
		for(int i=1; i<4; i++)
		{
			for(int a=0; a<4; a++)
			{
				e[i][a] = v[i][a];
			}

			for(int j=0; j<i; j++)
			{
				double dotprod = 0;
				double enorm = 0;

				// evaluate the dot products
				for(int a=0; a<4; a++)
					for(int b=0; b<4; b++)
					{
						dotprod += g[a][b]*v[i][a]*e[j][b];
						enorm += g[a][b]*e[j][a]*e[j][b];
					}
				// subtract the projections of the existing u's from this one
				for(int a=0; a<4; a++)
				{
					e[i][a] -= (dotprod/enorm) * e[j][a];
				}
			}
		}

		// flip signs as necessary to make a right-handed set
		if(e[1][1] < 0) for(int a=0; a<4; a++) e[1][a] *= -1;
		if(e[2][2] > 0) for(int a=0; a<4; a++) e[2][a] *= -1;
		if(e[3][3] < 0) for(int a=0; a<4; a++) e[3][a] *= -1;

		// swap vectors 1 and 3
		for(int a=0; a<4; a++)
		{
			vec[0][a] = e[0][a];
			vec[1][a] = e[3][a];
			vec[2][a] = e[2][a];
			vec[3][a] = e[1][a];
		}
	}

	void normalise()
	{
		// normalise the vectors
		for(int i=0; i<4; i++)
		{
			double enorm = 0;

			// evaluate the dot product
			for(int a=0; a<4; a++)
				for(int b=0; b<4; b++)
					enorm += g[a][b]*BasisVectors<T>::vectors[i][a]*BasisVectors<T>::vectors[i][b];

			for(int a=0; a<4; a++)
				BasisVectors<T>::vectors[i][a] /= sqrt(abs(enorm));
		}
	}
};

#endif //CUDAKERR_GRAMSCHMIDT_BASIS_H
