/** \file ImplicitUpwind.h
\brief  ImpilicitUpwind scheme declaration
*/

#pragma once
#include <vector> // std vector upon which our Vector is based
#include "matrix.h"
#include "vector.h"
#include "AbstractScheme.h"
///Class ImplicitUpwind for solving PDE, has access to public class called AbstractScheme
///ImplicitUpwind scheme is solved for BTCS( backward time, central space)
///This problem will be solved using LU factorization 
class ImplicitUpwind : public AbstractScheme
{
	Matrix * l;
	Matrix * u;
	Matrix * p;
	Matrix * prepareMatrixA(double cfl, double n);
public:
	ImplicitUpwind(double low, double up, double U, double N, double dx, double dt);
	~ImplicitUpwind();


	// Inherited via AbstractScheme
	virtual void calculateNewArray(double * uold, double * unew) override;

};

