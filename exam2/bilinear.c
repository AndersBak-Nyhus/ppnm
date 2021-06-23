#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include "binsearch.h"
//Interpolation in one direction
double linearInterpolation(double x1, double f_x1, double x2, double f_x2, double x)
{
	double result = (x - x1)/(x2-x1)*f_x2  + (x2-x)/(x2-x1)*f_x1;
	return result;
}


double bilinear(int nx, int ny, double* XData, double* YData, gsl_matrix* FData, double px, double py){
	//check which element in the vector is the closest with binsearch
	int ix = binsearch(nx, XData, px);

	int iy = binsearch(ny, YData, py);

	//get matrix elements corresponding to the nearest vector elements
	double F11 = gsl_matrix_get(FData,iy,ix);
	double F12 = gsl_matrix_get(FData,iy+1,ix);
	double F21 = gsl_matrix_get(FData,iy,ix+1);
	double F22 = gsl_matrix_get(FData,iy+1,ix+1);

	//interpolate along x
	double R1 =  linearInterpolation(XData[ix], F11, XData[ix+1],F21, px);
	double R2 =  linearInterpolation(XData[ix], F12, XData[ix+1], F22, px);
	//interpolate along y yields bilinear interpolation
	double  P =  linearInterpolation(YData[iy],  R1, YData[iy+1],  R2, py);
	return P;
}

