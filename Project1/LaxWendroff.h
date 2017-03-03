/** \file LaxWendroff.h
\brief  LaxWendroff scheme seclaration
*/
#pragma once
#include <vector> // std vector upon which our Vector is based
#include "AbstractScheme.h"

///Class LaxWendroff for solving PDE
class LaxWendroff : public AbstractScheme
{
public:
	LaxWendroff(double low, double up, double U, double N, double dx, double dt);
	~LaxWendroff();


	// Inherited via AbstractScheme
	virtual void calculateNewArray(double * uold, double * unew) override;

};

