/** \filRichtmyer scheme seclaration
*/
#pragma once
#include <vector> // std vector upon which our Vector is based
#include "AbstractScheme.h"

///Class Richtmyer for solving PDE
class RichtmyerMultistep : public AbstractScheme
{
public:
	RichtmyerMultistep(double low, double up, double U, double N, double dx, double dt);
	~RichtmyerMultistep();

	/**
	*Normal public method that returns a double
	*@return vectro policzHalf 
	*/
	double policzStepHalf(double*uold, double c, int i);



	// Inherited via AbstractScheme
	virtual void calculateNewArray(double * uold, double * unew) override;

};


