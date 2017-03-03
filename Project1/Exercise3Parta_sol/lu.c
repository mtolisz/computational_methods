
/* LU factorisation routine
* Takes in a matrix a of size n and produces the lower (l) and
* upper (u) triangular matrices that factorise a 
*/
void lu_fact(double** a, double** l, double** u, int n) 
{
	double** temp;
	double mult;
	int i,j,k;

	temp = matrix(n);

	// copy a to temp
	copy_matrix(a,temp,n);
	
	// LU (Doolittle's) decomposition without pivoting
	for (k = 0; k < n - 1; k++) {
		for (i = k + 1; i < n; i++) {
			if (fabs(temp[k][k]) < 1.e-07) 
			{ 
				printf("pivot is zero\n"); 
				exit(1);
			}
			mult = temp[i][k]/temp[k][k];
			temp[i][k] = mult;                      // entries of L are saved in temp
			for (j = k + 1; j < n; j++) { 
				temp[i][j] -= mult*temp[k][j];      // entries of U are saved in temp
				if (fabs(temp[i][i]) < 1.e-07)
				{ 
					printf("pivot is zero\n"); 
					exit(1);
				}
			}
		}
	}

	// create l and u from temp
	for (i=0; i<n; i++) l[i][i] = 1.0;
	for (i=1; i<n; i++) 
	  for (j=0; j<i; j++) l[i][j] = temp[i][j];

	for (i=0; i<n; i++)
		for (j=i; j<n; j++) u[i][j] = temp[i][j];

	free_matrix(temp,n);
} 


/*
* Solves the equation LUx = b by performing forward and backward
* substitution. Output is the solution vector x
*/
void lu_solve(double** l, double** u, double* b, int n, double* x)
{
	double* temp;
	int i,j; 

	temp = vector(n);

	// copy b to temp
	copy_vector(b,temp,n);
	
	// forward substitution for L y = b.
	for (i = 1; i < n; i++) 
		for (j = 0; j < i; j++) 
			temp[i] -= l[i][j]*temp[j];
	
  
	// back substitution for U x = y.  
	for (i = n - 1; i >= 0; i--) {
		for (j = i + 1; j < n; j++) temp[i] -= u[i][j]*temp[j];
		temp[i] /= u[i][i];
	}

	// copy solution into x
	for (i=0; i<n; i++) x[i] = temp[i];

	free_vector(temp,n);
}

/*
* Computes the permutation matrix P such that the matrix PA can be
* factorised into LU and the system PA = Pb can be solved by forward and 
* backward substitution. Output is the permutation matrix P
*/
void reorder(double** a, int n, double** p) 
{
// Note: pivoting information is stored in temperary vector pvt

	int i,j,k;
	int* pvt;
	int pvtk,pvti;
	double* scale;
	double aet, tmp, mult;
	double** temp;

	pvt = ivector(n);
	temp = matrix(n);

	// copy a into temp
	copy_matrix(a,temp,n);

	for (k = 0; k < n; k++) pvt[k] = k;

	scale = vector(n);             // find scale vector
	for (k = 0; k < n; k++) {
		scale[k] = 0;
		for (j = 0; j < n; j++) 
			if (fabs(scale[k]) < fabs(temp[k][j])) scale[k] = fabs(temp[k][j]);
	} 

	for (k = 0; k < n - 1; k++) {            // main elimination loop

	// find the pivot in column k in rows pvt[k], pvt[k+1], ..., pvt[n-1]
		int pc = k; 
		aet = fabs(temp[pvt[k]][k]/scale[k]);
		for (i = k + 1; i < n; i++) {
			tmp = fabs(temp[pvt[i]][k]/scale[pvt[i]]); 
			if (tmp > aet) {
				aet = tmp; 
				pc = i;
			}
		}
		if (fabs(aet) < 1.e-07) 
		{ 
			printf("pivot is zero\n"); 
			exit(1); 
		}
		if (pc != k) {                      // swap pvt[k] and pvt[pc]
			int ii = pvt[k];
			pvt[k] = pvt[pc];
			pvt[pc] = ii;
		}

		// now eliminate the column entries logically below mx[pvt[k]][k]
		pvtk = pvt[k];                           // pivot row
		for (i = k + 1; i < n; i++) {
			pvti = pvt[i];
			if (temp[pvti][k] != 0) {
				mult = temp[pvti][k]/temp[pvtk][k]; 
				temp[pvti][k] = mult;
				for (j = k + 1; j < n; j++) temp[pvti][j] -= mult*temp[pvtk][j];
			}
		}
	}
	for (i=0; i<n; i++) p[i][pvt[i]]=1.0;

	free_ivector(pvt,n);
	free_matrix(temp,n);
} 
