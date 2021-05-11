#include<assert.h>
double linterp(int n, double* x, double* y, double z){
	assert(n>1 && z>=x[0] && z<=x[n-1]);
	int i=0, j=n-1;
	while(j-i>1){
		int m=(i+j)/2;
		if(z>=x[m]) i=m;
		else j=m;
	}
	double ai=y[i];
	double bi=(y[i+1]-y[i])/(x[i+1]-x[i]);
	return ai+bi*(z-x[i]);
}

double linterp_integ(int n, double* x, double* y, double z){
        assert(n>1 && z>=x[0] && z<=x[n-1]);
        int i=0, j=n-1;
        while(j-i>1){
                int m=(i+j)/2;
                if(z>x[m]) i=m;
                else j=m;
        }
double Integral = 0;
double ak;
double bk;
for( int k=0; k <= i;k++){

	ak = y[k];
	bk = (y[k+1]-y[k])/(x[k+1]-x[k]);
	Integral += ak*(x[k+1]-x[k])+bk*(x[k+1]-x[k])*(x[k+1]-x[k])/2;//ak*(x[k+1]-x[k])+bk*(x[k+1]-x[k])*(x[k+1]-x[k])/2;
	
}
return Integral;
}
