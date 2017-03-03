#include "matrix.h"

void lu_fact(const Matrix& a, Matrix& l, Matrix& u, int n);
void lu_solve(const Matrix& l, const Matrix& u, const Vector& b, int n, Vector& x);
void reorder(const Matrix& a, int n, Matrix& p);

int main()
{
	int n;
	try {
		// input the size of the system
		std::cout << "input the size of the system\n";
		std::cin >> n;

		// declare the matrix a, vectors b and x
		Matrix a(n,n);
		Vector x(n), b(n);

		// input matrix a and right hand side vector b
		std::cin >> a;
		std::cin >> b;

		// declare l and u
		Matrix l(n,n), u(n,n);

		// declare p
		Matrix p(n,n);

		// reorder
		reorder(a,n,p);

		// multiply p by a
		Matrix m = p*a;

		// factorise
		lu_fact(m,l,u,n);

		// output l and u
		std::cout << l;
		std::cout << u;
	
		// solve
		lu_solve(l,u,p*b,n,x);

		// check solution

		if (a*x == b) {
			std::cout << "solution correct\n";
			// print solution
			std::cout << x;
		} else std::cout << "solution false\n";
	} catch (std::exception& e){
		// Catching other errors
		std::cerr << "std::exception caught" << std::endl;
		std::cerr << "Type: " << typeid(e).name() << std::endl;
		std::cerr << "What: " << e.what() << std::endl;
	}	
	return 0;
}