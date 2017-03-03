/** \file ExplicitUpwind.h
	\brief  ExplicitUpwind scheme seclaration
*/


#pragma once
#include <vector> // std vector upon which our Vector is based
#include "AbstractScheme.h"
///Class ExplicitUpwind for solving PDE
class ExplicitUpwind : public AbstractScheme
{
public:
	///CONSTRUCTOR
	///@param values type double for simulation 
	ExplicitUpwind(double low, double up, double U, double N, double dx, double dt);
	//DESTRUCTOR
	~ExplicitUpwind();

	///Calculate values for next step and stored into array
	/// Inherited via AbstractScheme
	virtual void calculateNewArray(double * uold, double * unew) override;

};

