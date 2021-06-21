#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "input2Array.h"
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include "bilinear.h"

int main(void)
{
	double numOfPts = 40;
	double *Xdata = malloc(numOfPts*sizeof(double));
	double *Ydata = malloc(numOfPts*sizeof(double));
	double numOfPts2 = 100;
	double *Fdata = malloc(numOfPts2*sizeof(double));
	double *Fdata1 = malloc(numOfPts2*sizeof(double));
	//read input.txt so as to create the vectors. collum 1 is x and collum 2 is y
	//the third line in the input.txt corresponds to the point P(x,y)
	// where we want to know the value of the function
	input2Array(Xdata,Ydata,"dataPts.txt");
	
	input2Array(Fdata,Fdata1,"dataVals.txt");

	int n = 5;
	gsl_matrix* F =gsl_matrix_alloc(n,n);
	gsl_matrix_set(F,0,0,Fdata[0]);
	gsl_matrix_set(F,1,0,Fdata[1]);
	gsl_matrix_set(F,2,0,Fdata[2]);
	gsl_matrix_set(F,3,0,Fdata[3]);
	gsl_matrix_set(F,4,0,Fdata[4]);
        gsl_matrix_set(F,0,1,Fdata[5]);
        gsl_matrix_set(F,1,1,Fdata[6]);
        gsl_matrix_set(F,2,1,Fdata[7]);
        gsl_matrix_set(F,3,1,Fdata[8]);
        gsl_matrix_set(F,4,1,Fdata[9]);
        gsl_matrix_set(F,0,2,Fdata[10]);
        gsl_matrix_set(F,1,2,Fdata[11]);
        gsl_matrix_set(F,2,2,Fdata[12]);
        gsl_matrix_set(F,3,2,Fdata[13]);
        gsl_matrix_set(F,4,2,Fdata[14]);
        gsl_matrix_set(F,0,3,Fdata[15]);
        gsl_matrix_set(F,1,3,Fdata[16]);
        gsl_matrix_set(F,2,3,Fdata[17]);
        gsl_matrix_set(F,3,3,Fdata[18]);
        gsl_matrix_set(F,4,3,Fdata[19]);
        gsl_matrix_set(F,0,4,Fdata[20]);
        gsl_matrix_set(F,1,4,Fdata[21]);
        gsl_matrix_set(F,2,4,Fdata[22]);
        gsl_matrix_set(F,3,4,Fdata[23]);
        gsl_matrix_set(F,4,4,Fdata[24]);
//printf("%g\n",Xdata[0]);
	double a = -1; double b =1;
//double c = 0.5;
	for(double px = a; px <= b; px += 0.2){
		for(double py = a; py <= b; py += 0.2){
			double fp = bilinear(5, Xdata, Ydata, F,px,py);
			printf("%g\t%g\t%g\n",px,py,fp);
		}
	
	}
	//Interpolation in the x-direction
//	double R1 =  linearInterpolation(Xdata[0], F11, Xdata[1], F21, Xdata[2]);
//	double R2 =  linearInterpolation(Xdata[0], F12, Xdata[1], F22, Xdata[2]);
	//by doing interpolation first in the x-direction on the two lines and then over y, with the two interpolated values, the total sequence is a bilinear interpolation rutine


	free(Xdata);
	free(Ydata);
	free(Fdata);
	free(Fdata1);
	gsl_matrix_free(F);
}
