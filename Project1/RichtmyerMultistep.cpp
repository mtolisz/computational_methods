#include "RichtmyerMultistep.h"


RichtmyerMultistep::RichtmyerMultistep(double low, double up, double U, double N, double dx, double dt)
	:AbstractScheme(low, up, U, N, dx, dt)
{
}

RichtmyerMultistep::~RichtmyerMultistep()
{
}

double RichtmyerMultistep::policzStepHalf(double * uold, double c, int i)
{
	return 0.5 * (uold[i + 1] + uold[i - 1]) - 0.25 * c * (uold[i + 1] - uold[i - 1]);
}

void RichtmyerMultistep::calculateNewArray(double * uold, double * unew)
{
	double c = this->U * this->dt / this->dx;
	for (int i = 2; i < N - 2; i++)
	{
		unew[i] = uold[i] - 0.5 * c * (policzStepHalf(uold, c, i + 1) - policzStepHalf(uold, c, i - 1));
	}

	unew[0] = uold[0];
	unew[(int)N - 1] = uold[(int)N - 1];

	unew[1] = uold[1];
	unew[(int)N - 2] = uold[(int)N - 2];
	for (int i = 0; i < N; i++)
	{
		uold[i] = unew[i];
	}
}


