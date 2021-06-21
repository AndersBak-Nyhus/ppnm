#ifndef HAVE_BILINEAR_H
#define HAVE_BILINEAR_H
#include<gsl/gsl_matrix.h>

double bilinear(int n, double* XData, double* YData, gsl_matrix* FData, double px, double py );

#endif 
