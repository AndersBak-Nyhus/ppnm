#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "ann.h"


//Function along with derivative and antiderivative
//Gaussian Wavelet
double activation_func(double x){
	return x*exp(-x*x);
}
double activation_der(double x){
	return (1-2*x*x)*exp(-x*x);
}
double activation_integ(double x){
	return -exp(-x*x)*1./2;
}

//Interpolation

double func(double x){
	return sin(x)*exp(-x);
}
//Derivative

double func_der(double x){
	return (cos(x)-sin(x))*exp(-x);
}
//Antiderivative

double func_integ(double x, double x0){
	return -1./2*(sin(x)+cos(x))*exp(-x)+1./2*(sin(x0)+cos(x0))*exp(-x0);
}

int main(){
	int nx = 200;
	gsl_vector* xs = gsl_vector_alloc(nx);
	gsl_vector* ys = gsl_vector_alloc(nx);
	gsl_vector* yms = gsl_vector_alloc(nx);
	gsl_vector* Ys = gsl_vector_alloc(nx);

	int n = 3; //number of neurons
	ann* network = ann_alloc(n,activation_func, activation_der, activation_integ);
	double xmin = 0;
	double xmax = 2; // interval on x-axis
	// set parameters
	for(int i=0; i<network->n; i++){
		gsl_vector_set(network -> params, 3*i,xmin + (xmax - xmin) * i/(network -> n-1));
		gsl_vector_set(network -> params, 3*i+1,1);
		gsl_vector_set(network -> params, 3*i+2,1);
	}
	for(int i = 0; i < nx; i++){
		gsl_vector_set(xs,i,xmin + (xmax - xmin)*i/(nx-1));
		gsl_vector_set(ys,i,func(gsl_vector_get(xs,i)));
		gsl_vector_set(yms,i,func_der(gsl_vector_get(xs,i)));
		gsl_vector_set(Ys,i,func_integ(gsl_vector_get(xs,i),xmin));
	}

	ann_train(network, nx, xs, ys);
	printf("Part A)\n");
	printf("Activation function: x*exp(-x^2)\n\n");
	for(int i = 0; i<network->n; i++){
		double ai = gsl_vector_get(network->params,3*i);
		double bi = gsl_vector_get(network->params,3*i+1);
		double wi = gsl_vector_get(network->params,3*i+2);
		printf("i=%i, ai, bi, wi = %g %g %g\n", i, ai, bi, wi);
	}
	printf("\n\nPart B)\n");
	printf("See plot.png\n");
	FILE* Pts =fopen("Pts.txt", "w");
	for(int i=0; i<nx;i++) fprintf(Pts, "%g %g %g %g\n", gsl_vector_get(xs,i), gsl_vector_get(ys,i), gsl_vector_get(yms,i), gsl_vector_get(Ys,i));

	// functions
	FILE* Funcs = fopen("Funcs.txt", "w");
	double dz=0.2;
	for(double z = xmin; z<=xmax; z+=dz){
		fprintf(Funcs, "%g %g %g %g\n", z, ann_response(network,z), ann_der(network,z), ann_integ(network,z,xmin));
	}
	ann_free(network);
	fclose(Pts);
	fclose(Funcs);


	gsl_vector_free(xs);
	gsl_vector_free(ys);
	gsl_vector_free(yms);
	gsl_vector_free(Ys);

	return 0;
}
