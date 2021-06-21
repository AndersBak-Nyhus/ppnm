#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "input2Array.h"

double linearInterpolation(double x1, double f_x1, double x2, double f_x2, double x)
{
	double result = (x - x1)/(x2-x1)*f_x2  + (x2-x)/(x2-x1)*f_x1;
	return result;
}

int main(void)
{
	double numOfPts = 10;
	double *Xdata = malloc(numOfPts*sizeof(double));
	double *Ydata = malloc(numOfPts*sizeof(double));

	//read input.txt so as to create the vectors. collum 1 is x and collum 2 is y
	//the third line in the input.txt corresponds to the point P(x,y)
	// where we want to know the value of the function
	input2Array(Xdata,Ydata,"x_i_y_i.txt");

	//function values at the four points first collum in inputVals.txt correspond to y1 and x1,x2 and the second collum y2 and x1,x2
	double *Fvalsy1 = malloc(numOfPts*sizeof(double));
	double *Fvalsy2 = malloc(numOfPts*sizeof(double));

	input2Array(Fvalsy1,Fvalsy2,"F_ijVals.txt");
	//double F11 =1.5; double F12 = 3.3; double F21 = 2.1; double F22 = 4.5;
	double F11 = Fvalsy1[0]; double F21 = Fvalsy1[1]; double F12 = Fvalsy2[0]; double F22 = Fvalsy2[1];

	//Interpolation in the x-direction
	double R1 =  linearInterpolation(Xdata[0], F11, Xdata[1], F21, Xdata[2]);
	double R2 =  linearInterpolation(Xdata[0], F12, Xdata[1], F22, Xdata[2]);
	printf("vector elements:\n");
	printf("x1=%g, x2=%g\n",Xdata[0],Xdata[1]);
	printf("y1=%g, y2=%g\n",Ydata[0],Ydata[1]);

	printf("\n point we want to interpolate:\n");
	printf("P=(x=%g,y=%g)\n",Xdata[2],Ydata[2]);

	printf("\n The interpolated value of the function at R1=(x,y1) is %g\n", R1);
	printf("The interpolated value of the function at R2=(x,y2) is %g\n", R2);

	//interpolation along y, giving the desired result at the point P
	double  P =  linearInterpolation(Ydata[0],  R1, Ydata[1],  R2, Ydata[2]);
	
	printf("\nThe interpolated value of the function at (x=%g, y=%g) is equal to %g\n\n", Xdata[3], Ydata[3], P);
	
	free(Xdata);
	free(Ydata);
	free(Fvalsy1);
	free(Fvalsy2);
	printf("\n change values of x_i,x_j and y_i,y_j in x_i_y_i.txt file as well as interpolated point in the third line\n");
	printf("\n change matrix values F_ij in F_ijVals.txt");
}
