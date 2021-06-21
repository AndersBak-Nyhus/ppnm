#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include "binsearch.h"

double linearInterpolation(double x1, double f_x1, double x2, double f_x2, double x)
{
	double result = (x - x1)/(x2-x1)*f_x2  + (x2-x)/(x2-x1)*f_x1;
	return result;
}


void bilinear(int n, double* XData, double* YData, gsl_matrix* FData, double px, double py){

	int ix = binsearch(n, XData, px);
	int k = n;
	int iy = binsearch(k, YData, py);
//	printf("integers: %d\t%d\n",ix,iy);
//	int ix = 1;
	double F11 = gsl_matrix_get(FData,ix,iy);
	double F12 = gsl_matrix_get(FData,ix,iy+1);
	double F21 = gsl_matrix_get(FData,ix+1,iy);
	double F22 = gsl_matrix_get(FData,ix+1,iy+1);

	double R1 =  linearInterpolation(XData[ix], F11, XData[ix+1],F21, px);
	double R2 =  linearInterpolation(XData[ix], F12, XData[ix+1], F22, px);

	double  P =  linearInterpolation(YData[iy],  R1, YData[iy+1],  R2, py);
	//printf("%g\t%g\t%g\n",px,py,P);
}

