#include "LaxWendroff.h"

LaxWendroff::LaxWendroff(double low, double up, double U, double N, double dx, double dt)
	: AbstractScheme(low, up, U, N, dx, dt)
{
}

LaxWendroff::~LaxWendroff()
{
}

void LaxWendroff::calculateNewArray(double * uold, double * unew)
{
	double c = this->U * this->dt / this->dx;
	for (int i = 1; i < N - 1; i++)
	
	unew[i] = uold[i] - 0.5 * c * (uold[i+1] - uold[i - 1]) + 0.5 * c * c * (uold[i + 1] - 2 * uold[i] + uold[i - 1]);
	unew[0] = uold[0];
	unew[(int)N - 1] = uold[(int)N - 1];

	for (int i = 0; i < N; i++)
	{
		uold[i] = unew[i];
	}
}
