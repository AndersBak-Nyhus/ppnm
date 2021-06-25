#include <stdio.h>
#include <math.h>
#include <complex.h>

complex plainmc(int dim, double f(double* x), double* a, double* b, int N);
complex quasimc(int dim, double f(double* x), double* a, double* b, int N);

double f(double * x){
	return 1./(M_PI*M_PI*M_PI) * 1./(1. - cos(x[0]) * cos(x[1]) * cos(x[2]));
}

double g(double* x){
	return 1./sqrt(x[0]);
}
double y(double* x){
	return 1./exp(1-x[0]);
}

double c(double* p){
	double x = p[0];
	double y = p[1];
	double R = 0.9;
	if(x*x + y*y < R*R) return 1;
	else return 0;
	}


int main(){
	double a[] = {0,0,0};
	double b[] = {M_PI, M_PI, M_PI};
	int N = 1e7;

	complex result = plainmc(3,f,a,b,N);

	//Compairing results
	printf("Part A)\n\n");
	printf ("Numerical with errors = %g +/- %g*i \n", creal(result), cimag(result));
	printf ("Analytical= %g \n", pow(tgamma(1./4), 4)/(4*M_PI*M_PI*M_PI));


	double a_test[] = {0};
	double b_test[] = {M_PI};
	complex result_test1b = plainmc(1,g,a_test,b_test,N);
	printf("\nTest functions: 1/sqrt(x) and 1/exp(1-x)\n\n");
	printf("1/sqrt(x) from 0 to pi = %g +/- %g*i\n", creal(result_test1b),cimag(result_test1b));


	complex result_test1c = plainmc(1,y,a_test,b_test,N);
	printf("1/exp(1-x) from 0 to pi = %g +/- %g*i\n",creal(result_test1c),cimag(result_test1c));


	//Part B)
	N = 1e6;
	printf("\nPart B)\n\n");
	int dim = 3;
	complex result2a = quasimc(dim, f, a, b, N);
	printf("Solution and error = %g +/- %g*i\n", creal(result2a), cimag(result2a));
	complex result2b = quasimc(1,g,a_test,b_test,N);
	printf("1/sqrt(x) from 0 to pi: %g +/- %g*i\n",creal(result2b), cimag(result2b));
	complex result2c = quasimc(1,y,a_test,b_test,N);
	printf("1/exp(1-x) from 0 to pi: %g +/- %g*i\n",creal(result2c),cimag(result2c));

	//Error scaling
	FILE* error = fopen("error.txt","w");
	double a_dim[] = {0,0};
	double b_dim[] = {2*M_PI, 2*M_PI};
	for(int N = 1000; N < 100000; N+= 1000){
		complex result_plain = plainmc(2,c,a_dim,b_dim,N);
		complex result_quasi = quasimc(2,c,a_dim,b_dim,N);
		fprintf(error, "%20d %20g %20g\n", N, cimag(result_plain), cimag(result_quasi));
	}
	printf("Scaling in error.png\n");

	fclose(error);






	return 0;
}
