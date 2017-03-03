#include "Initialization.h"
#include <cmath> // File math operations
double Initialization::sign(double x)
{
	if (x < 0) return -1;
	if (x > 0) return 1;
	return 0;
}

void Initialization::initSign(double *arr, double dX) {
	double x = low;

	for (int i = 0; i < N; i++, x += dX) {
		arr[i] = 0.5 * (sign(x) + 1.0);
	}
}

void Initialization::initExp(double * arr, double dX) {

		double x = low;
		for (int i = 0; i < N; i++, x += dX) {
			arr[i] = 0.5 * exp(-x*x);
		}

}






	






