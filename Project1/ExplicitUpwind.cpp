/*! \file ExplicitUpwind.cpp
\brief  Solution for the Explicit scheme for advection eqations*/
#include "ExplicitUpwind.h"

//!Constructor
ExplicitUpwind::ExplicitUpwind(double low, double up, double U, double N, double dx, double dt)
	: AbstractScheme(low, up, U, N, dx, dt)
{
}
//! Destructor
ExplicitUpwind::~ExplicitUpwind()
{
}
//! CalculateNewArray method of ExplicitUpwind
///return values of next step
void ExplicitUpwind::calculateNewArray(double * uold, double * unew)
{
	double c = this->U * this->dt / this->dx;
	for (int i = 1; i < N - 1; i++)
		unew[i] = uold[i] - c*(uold[i] - uold[i - 1]);
	unew[0] = uold[0];
	unew[(int)N-1] = uold[(int)N-1];
	for (int i = 0; i < N; i++)
	{
		uold[i] = unew[i];
	}
}
