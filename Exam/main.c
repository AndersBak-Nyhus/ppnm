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
	//numofPts is the number of points the function input2Array is expecting, but it doesn't always work on the exact number of pts, so in this case they are raised above
	double numOfPts = 40;
	//allocate memory for the x and y data
	double *Xdata = malloc(numOfPts*sizeof(double));
	double *Ydata = malloc(numOfPts*sizeof(double));

	double numOfPts2 = 100;
	//allocate memory for the function values
	double *Fdata = malloc(numOfPts2*sizeof(double));
	double *Fdata1 = malloc(numOfPts2*sizeof(double));
	//read dataPts.txt so as to create the vectors. collum 1 is x and collum 2 is y
	//function values are in dataVals.txt where the first collum corresponds to the point P(x_i,y_j)
	// where we want to know the value of the function
	input2Array(Xdata,Ydata,"dataPts.txt");
	
	input2Array(Fdata,Fdata1,"dataVals.txt");
	//The -3 comes from the program i used to make the txt files, which creates "invisible" data
	int nx = sizeof(Xdata)-3;
	int ny = sizeof(Ydata)-3;

	//create matrix
	gsl_matrix* F =gsl_matrix_alloc(nx,ny);
	int k = 0;
	for(int i = 0; i < nx; i += 1){
		for(int j = 0; j < ny; j +=1){
			int k = k+1; 
			gsl_matrix_set(F,j,i,Fdata[k-1]);
		}
	}



	//numbers needed in the sum below where a is the minimum in the data set and b is max
	double a = -1; double b =1;
	//sum over both x and y in steps of 0.1
	for(double px = a; px <= b; px += 0.1){
		for(double py = a; py <= b; py += 0.1){
			//send x and y to bilinear function to interpolate
			double fp = bilinear(nx, ny, Xdata, Ydata, F,px,py);
			//print data in output.txt to be able to plot it
			printf("\n%g\t%g\t%g\n",px,py,fp);
		}
	
	}
	//free allocated memory

	free(Xdata);
	free(Ydata);
	free(Fdata);
	free(Fdata1);
	gsl_matrix_free(F);
}
