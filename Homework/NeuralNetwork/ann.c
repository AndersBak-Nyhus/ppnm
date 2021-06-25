#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include "minimization.h"
#include "ann.h"

ann* ann_alloc(int n, double(*f)(double), double(*fm)(double), double(*F)(double)){ // function to allocate memory of network
	ann* network = malloc(sizeof(ann)); // initialize with zeros
	network -> n = n;
	network -> f = f;
	network -> fm = fm;
	network -> F=F;
	network -> params = gsl_vector_alloc(3*n); // three parameters (a,b,w) for each neuron
	return network;
}

//free allocated memory
void ann_free(ann* network){
	free(network -> params);
	free(network);
}


// Response function (Fp(x) in notes)
double ann_response(ann* network, double x){
	double sum = 0;
	for(int i = 0; i < network -> n ;i++){
		double a = gsl_vector_get(network -> params, 3*i+0);
		double b = gsl_vector_get(network -> params, 3*i+1);
		double w = gsl_vector_get(network -> params, 3*i+2);
		sum += network -> f((x-a)/b)*w;
	}
	return sum;
}


//derivative of ann function
double ann_der(ann* network, double x){
	double sum = 0;
	for(int i = 0; i < network -> n; i++){
		double a = gsl_vector_get(network -> params,3*i+0);
		double b = gsl_vector_get(network -> params,3*i+1);
		double w = gsl_vector_get(network -> params,3*i+2);
		sum += network -> fm((x-a)/b)*w/b;
	}
	return sum;
}

//The anti-derivative
double ann_integ(ann* network, double x, double x0){
	double sum = 0;
	for(int i = 0; i < network -> n; i++){
		double a = gsl_vector_get(network -> params,3*i);
		double b = gsl_vector_get(network -> params,3*i+1);
		double w = gsl_vector_get(network -> params,3*i+2);
		sum += network -> F((x-a)/b)*w*b - network -> F((x0-a)/b)*w*b;
	}
	return sum;
}




static int N;
gsl_vector* X;
gsl_vector* Y;


static ann* network_p;


double cost_func(gsl_vector* p){

	gsl_vector_memcpy(network_p -> params, p);

	double sum = 0;
	for(int k = 0; k < N; k++){
		double xk = gsl_vector_get(X,k);
		double yk = gsl_vector_get(Y,k);

		double fk = ann_response(network_p, xk);
		sum += (fk - yk)*(fk - yk);
	}

	return sum/N;

}
void ann_train(ann* network, int nx, gsl_vector* xs, gsl_vector* ys){
	network_p = network;
	X = gsl_vector_alloc(xs -> size);
	Y = gsl_vector_alloc(ys -> size);
	gsl_vector* p = gsl_vector_alloc(network->params -> size);
	//printf("%ld, %ld, %ld, %ld, %ld\n",X->size, Y->size, xs->size, ys->size,p->size);
	N = nx;
	// X = xs, Y = ys
	for (int i = 0; i < xs -> size; i++){
		gsl_vector_set(X, i, gsl_vector_get(xs,i));
		gsl_vector_set(Y,i,gsl_vector_get(ys,i));
	}
	//int d = 3*network -> n;
	double acc = 1e-3;
	gsl_vector_memcpy(p,network->params);

	quasi_newtonMethod(cost_func, p, acc);
	gsl_vector_memcpy(network->params,p);

	gsl_vector_free(p);
}




