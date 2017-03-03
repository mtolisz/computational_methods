#pragma once
#include <vector> // std vector upon which our Vector is based
#include "matrix.h" ////we use matrix in Matrix code
#include "vector.h" //we use vector in Matrix code

/**
*  A AbstractScheme class for data storage of a 2D array of doubles
*  \n The implementation is derived from the standard container vector std::vector
*	\n The implementation is derived from the standard container vector std::vector
*  \n We use private inheritance to base our vector upon the library version whilst
*  \nallowing usto expose only those base class functions we wish to use - in this
*  \ncase the array access operator []
*
* The Matrix class provides:
* \n-basic constructors for creating a matrix object from other matrix object,
* \nor by creating empty matrix of a given size,
* \n-input and oput operation via >> and << operators using keyboard or file
* \n-basic operations like access via [] operator, assignment and comparision
*/

class AbstractScheme
{
protected:
	double low; //lower boundary
	double up; //upper boundary
	double U; // condition 
	double N; //grid size
	double dx; // space step
	double dt; //time step 

public:
	//!CONSTRUCTOR
	/*!
	values are types of double 
	*/

	explicit AbstractScheme(double low, double up, double U, double N, double dx, double dt);
	//! DESTRUCTOR
	virtual ~AbstractScheme();
	//! Constructor of new method 
	virtual void calculateNewArray(double*uold, double*unew) = 0;
};

