#include "ImplicitUpwind.h"
#include "matrix.h"
#include "vector.h"

void lu_fact(const Matrix& a, Matrix& l, Matrix& u, int n);
void lu_solve(const Matrix& l, const Matrix& u, const Vector& b, int n, Vector& x);
void reorder(const Matrix& a, int n, Matrix& p);

Matrix * ImplicitUpwind::prepareMatrixA(double cfl, double n)
{
	Matrix * a = new Matrix(n, n);
	double alfa = 0.5 * cfl;

	unsigned int index = 0;

	(*a)[0][0] = 1.0;
	(*a)[0][1] = alfa;
	for (unsigned int row = 1; row < n - 1; row++)
	{
		(*a)[row][index] = -alfa;
		(*a)[row][index + 1] = 1.0;
		(*a)[row][index + 2] = alfa;
		index++;

	}
	(*a)[n - 1][n - 2] = -alfa;
	(*a)[n - 1][n - 1] = 1.0;

	return a;
}

ImplicitUpwind::ImplicitUpwind(double low, double up, double U, double N, double dx, double dt)
	: AbstractScheme(low, up, U, N, dx, dt)
{
	double n = N - 2;
	double cfl = U * dt / dx;
	Matrix * a = prepareMatrixA(cfl, n);

	// declare p
	this->p = new Matrix(n, n);
	this->l = new Matrix(n, n);
	this->u = new Matrix(n, n);
	// reorder
	reorder(*a, n, *p);

	// multiply p 0by a00
	Matrix m = (*p)*(*a);

	// factorise
	lu_fact(m, (*l), (*u), n);
	
	delete a;
}

ImplicitUpwind::~ImplicitUpwind()
{
	delete l;
	delete u;
	delete p;
}

void ImplicitUpwind::calculateNewArray(double * uold, double * unew)
{
	double n = N - 2;
	double c = this->U * this->dt / this->dx;
	Vector x(n), b(n);

	for (size_t i = 1; i < N - 1; i++)
	{
		b[i - 1] = uold[i];
	}

	//poprawka
	b[0] += 0.5 * c * uold[0];
	b[n - 1] -= 0.5 * c * uold[(int)N - 1];
	//koniec_poprawki

	lu_solve((*l), (*u), (*p)*b, n, x);

	for (size_t i = 1; i < N  - 1; i++)
	{
		uold[i] = x[(i - 1)];
	}
}
