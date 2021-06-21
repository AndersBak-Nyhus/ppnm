#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double linearInterpolation(double x1, double f_x1, double x2, double f_x2, double x)
{
	double result = (x - x1)/(x2-x1)*f_x2  + (x2-x)/(x2-x1)*f_x1;
	return result;
}


void bilinear(int n, double* XData, double* YData, double* FData, double px, double py){

	int ix = binsearch(n, XData, px);
	int iy = binsearch(n, YData, py);



	double R1 =  linearInterpolation(XData[ix], FData[ix,ix], XData[ix+1], FData[ix+1,ix], XData[ix]+px);
	double R2 =  linearInterpolation(XData[ix], FData[ix,ix+1], XData[ix+1], FData[ix+1,ix+1], XData[ix]+px);

	double  P =  linearInterpolation(YData[iy],  R1, YData[iy+1],  R2, YData[iy]+py);
	return P;
}

