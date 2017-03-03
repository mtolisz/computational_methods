#pragma once
/**	\ file	Initialization.h
*	\ brief	Solving  Advection equation analytically
*	Signum Analytical and Exponential functions
*	These two solutions are  found  here 
*/

class Initialization
{
private:

public:
	const double u = 1.75;
	const double low = -50.0; //low boundary
	const double up = +50.0; //upper boundary
	int N = 100; //number of points

/**	Called method initSign 
 *	@param	*arr	pointer to Array.
 *	@param	dX	Point in space.
 */
	void initSign(double * arr, double dX);

	/**	Called method initExp 
*	@param	*arr	pointer to Array.
*	@param	dX	Point in space.
*/
	void initExp(double * arr, double dX);

	/**	Analytical solution  for	eqaution  with	sign  function .
	*	@param	value (double)
	*	@param	N (integer)
	*/
	double sign(double value);
	Initialization(int N)
	{
		this->N = N;
	}


};