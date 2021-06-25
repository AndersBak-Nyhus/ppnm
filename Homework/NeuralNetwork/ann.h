#include "gsl/gsl_vector.h"
#ifndef HAVE_ANN_H
#define HAVE_ANN_H
typedef struct{
	int n;
	double(*f)(double);
	double(*fm)(double);
	double(*F)(double);
	gsl_vector* params;

} ann;

//functions
ann* ann_alloc(int n, double (*f)(double), double (*fm)(double), double (*F)(double));
void ann_free(ann* network);
double ann_response(ann* network, double x);
double ann_der(ann* network, double x);
double ann_integ(ann* network, double x, double x0);
double cost_func(gsl_vector* p);
void ann_train(ann* network, int nx, gsl_vector* xs, gsl_vector* ys);

#endif
