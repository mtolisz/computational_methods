#include <iostream> // std::cout
#include <iterator> // std::iterator
#include <fstream> //File IO operations
#include <cmath> // File math operations
#include "Initialization.h"
#include "ExplicitUpwind.h"
#include "ImplicitUpwind.h"   
#include "LaxWendroff.h"
#include "AbstractScheme.h"
#include "RichtmyerMultistep.h"

const double U = 1.75;
const double LOW = -50; //low boundary
const double UP = 50; //low boundary

void runScheme(AbstractScheme * scheme, double u, double low, double up, int n, double dt, double dx, const char * filename, bool isSign);
double loadDT();
int loadn();
int loading();
void checkCFL(double cfl, int option);
void sayDT(double dx, double a, int option);
int main() {

	int option = loading();

	int n = loadn();
	Initialization initialization(n);

	double dx = (UP - LOW) / n;
	sayDT(dx, U, option);
	double dt = loadDT();
	std::cout << "delta(x) = " << dx << std::endl;

	AbstractScheme * scheme = NULL;

	checkCFL(U*dt / dx, option);

	switch (option)
	{
	case 1:
		scheme = new ExplicitUpwind(LOW, UP, U, n, dx, dt);
		runScheme(scheme, U, LOW, UP, n, dt, dx, "EXPresultExp.csv", false);
		runScheme(scheme, U, LOW, UP, n, dt, dx, "EXPresultSign.csv", true);
		break;
	case 2:
		scheme = new ImplicitUpwind(LOW, UP, U, n, dx, dt);
		runScheme(scheme, U, LOW, UP, n, dt, dx, "IMPresultExp.csv", false);
		runScheme(scheme, U, LOW, UP, n, dt, dx, "IMPresultSign.csv", true);
		break;
	case 3:
		scheme = new LaxWendroff(LOW, UP, U, n, dx, dt);
		runScheme(scheme, U, LOW, UP, n, dt, dx, "LAXWEnresultExp.csv", false);
		runScheme(scheme, U, LOW, UP, n, dt, dx, "LAXWEnresultSign.csv", true);
		break;
	case 4:
		scheme = new RichtmyerMultistep(LOW, UP, U, n, dx, dt);
		runScheme(scheme, U, LOW, UP, n, dt, dx, "RICHTresultExp.csv", false);
		runScheme(scheme, U, LOW, UP, n, dt, dx, "RICHTresultSign.csv", true);
		break; 
	case 0: std::cout << "Thank you"; break;
	default: std::cout << "Unknown option";
	}
	if(scheme !=NULL)
		delete scheme;
//	system("pause");
	return 0;
}

void runScheme(AbstractScheme * scheme, double u, double low, double up, int n, double dt, double dx, const char * filename, bool isSign)
{
	Initialization initialization(n);

	double * unew = new double[n];
	double * anal = new double[n];
	double * uold = new double[n];
	double * analytical = new double[n];
	double c_upwind = (u * dt / dx);
	int timeStops[] = { 5, 10, 20 };
	double time = 0.0;


	//Initial values for T = 0
	if (isSign)
		initialization.initSign(uold, dx);
	else
		initialization.initExp(uold, dx);

	std::ofstream result;
	result.open(filename);

	//nORMS
	for (int t = 0; t < 3; t++)
	{
		while (time < timeStops[t])
		{
			time += dt;
			scheme->calculateNewArray(uold, unew);
		}
		double x = low, normSum = 0.0;
		for (int i = 0; i < n; i++, x += dx)
		{
			double a;
			if (isSign)
				a = 0.5 * (initialization.sign(x - u* time) + 1.0);
			else
				a = 0.5 * exp(-(x - u * time) * (x - u * time));

			anal[i] = a;
			normSum += (a - uold[i]) * (a - uold[i]);
		}
		std::cout << "T = " << timeStops[t] << " (" << time << ")" << std::endl;
		std::cout << "norm 2 = " << sqrt(normSum) << std::endl;

		for (int i = 0; i < n; i++) {
			result << uold[i] << ",";
		}
		result << std::endl;

		for (int i = 0; i < n; i++) {
			result << anal[i] << ",";
		}
		result << std::endl;
	}

	result.close();

	delete[] unew;
	delete[] uold;
	delete[] anal;
}

int loading()
{
	int nOption;
	std::cout << "Choose scheme" << std::endl;
	std::cout << "1. Explicit Upwind" << std::endl;
	std::cout << "2. Implicit Upwind" << std::endl;
	std::cout << "3. Lax - Wendroff" << std::endl;
	std::cout << "4. Richtmyer multi-step" << std::endl;     
	std::cout << "0. Exit" << std::endl;
	std::cout << "Your choose: ";
	std::cin >> nOption;



	switch (nOption)
	{
	case 1: std::cout << "ExplicitUpwind" << std::endl; break;//<< Explicit.OUTPUT(f0) << std::endl;
	case 2: std::cout << "Implicit Upwind" << std::endl; break;//<< Implicit.OUTPUT(f0) << std::endl;
	case 3: std::cout << "Lax - Wendroff" << std::endl; break;//<< LaxW.OUTPUT(f0) << std::endl;
	case 4:std::cout << "Richtmyer multi-step" << std::endl; break; //<< Richtmyer.OUTPUT(f0) << std::endl;
	case 0: std::cout << "Thank you" << std::endl; break;
	default: std::cout << "Unknown option!" << std::endl;
	}

	return nOption;
}

void checkCFL(double cfl, int option)
{
	switch (option)
	{
	case 1:
		if (cfl >= 1 || cfl <= 0)
		{
			std::cout << "Incorrect CFL, schema becomes unstable!" << std::endl;
			exit(-1);
		}
		break;
	case 2:
		if (cfl <= 0)
		{
			std::cout << "Incorrect CFL, schema becomes unstable!" << std::endl;
			exit(-1);
		}
		break;
	case 3:
		if (cfl >= 1 || cfl <= 0)
		{
			std::cout << "Incorrect CFL, schema becomes unstable!" << std::endl;
			exit(-1);
		}
		break;
	case 4:
		if (cfl >= 2 || cfl <= 0)
		{
			std::cout << "Incorrect CFL, schema becomes unstable!" << std::endl;
			exit(-1);
		}
		break;
	default: 
		exit(-1);
	}
}

void sayDT(double dx, double a, int option)
{
	switch (option)
	{
	case 1:

		std::cout << "dt should be <= " << dx / a<< std::endl;
		
		break;
	case 2:
		std::cout << "dt should be >=0 " << std::endl;

		break;
	case 3:
		std::cout << "dt should be <= " << dx / a << std::endl;

		break;
	case 4:
		std::cout << "dt should be <= " << 2* dx / a << std::endl;

		break;
	default:
		exit(-1);
	}
}

double loadDT()
{
	double dt;
	std::cout << "Enter dt\n";
	std::cin >> dt;
	return dt;
}

int loadn()
{
	int n;
	std::cout << "Enter grid points\n";
	std::cin >> n;
	return n;
}