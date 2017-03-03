/**\File AbstractScheme.h
*\ This source file contains on implementation of the AbstractScheme class
*/
#include "AbstractScheme.h"
//! CONSTRUCTOR
AbstractScheme::AbstractScheme(double low, double up, double U, double N, double dx, double dt)
	:low(low), up(up), N(N), dx(dx), U(U), dt(dt)
{
}
//! DESTRUCTOR
AbstractScheme::~AbstractScheme()
{
}
