#include<stdio.h>
#include<stdlib.h>
#define RND (double)rand()/RAND_MAX
double linterp(int,double*,double*,double);
#include<math.h>

int main(){
	int n=10;
	double x[n],y[n];
	for(int i=0;i<n;i++){
		x[i] = 0.5+i;
		y[i] = 0.2*i+RND;
		printf("%g %g\n",x[i],y[i]);
	}
	printf("\n\n");
	double dz = 0.01;
	for(double z=x[0];z<=x[n-1];z+=dz){
		double lz = linterp(n,x,y,z);
		printf("%g %g\n",z,lz);
	}
	double a=-9.9,b=9.9;
	n/=2;
	fprintf(stderr,"#m=0,S=4\n");
	for(int i=0;i<n;i++){
		x[i]=a+(b-a)*i/(n-1);
		y[i]=x[i]*x[i]*x[i];
		fprintf(stderr,"%g %g\n",x[i],y[i]);
	}
	fprintf(stderr,"#m=1,S=0\n");
	for(double z=x[0];z<=x[n-1];z+=dz) fprintf(stderr,"%g %g\n",z,linterp(n,x,y,z));

return 0;
}
