#include<stdio.h>
#include<stdlib.h>
#define RND (double)rand()/RAND_MAX
double linterp(int,double*,double*,double);
#include<math.h>
double linterp_integ(int,double*,double*,double);
#include<gsl/gsl_interp.h>


int main(int argc,char* argv[]){
FILE* outFileStream_lin = fopen(argv[1],"w");
//creating random points to interpolate
	int n=20;
	double x[n],y[n];
	for(int i=0;i<n;i++){
		x[i] = 0.5+i;
		y[i] = 0.2*i+RND;
		printf("%g %g\n",x[i],y[i]);
	}
	printf("\n\n\n");

	gsl_interp* linear = gsl_interp_alloc(gsl_interp_linear,n); 

	gsl_interp_init(linear,x,y,n);

	double dz = 0.01;
	for(double z=x[0];z<=x[n-1];z+=dz){
		double lz = linterp(n,x,y,z);
		
		double integral_lz = linterp_integ(n,x,y,z);
		
		double gsl_interp_l = gsl_interp_eval(linear,x,y,z,NULL);
		
		double gsl_integ_l = gsl_interp_eval_integ(linear,x,y,x[0],z,NULL);
		
		printf("%g %g %g %g %g\n",z,lz,integral_lz,gsl_interp_l,gsl_integ_l);
	}
//	double a=-9.9,b=9.9;
	n/=2;
	fprintf(stderr,"#m=0,S=4\n");
//	for(int i=0;i<n;i++){
//		x[i]=a+(b-a)*i/(n-1);
//		y[i]=x[i]*x[i]*x[i];
//		fprintf(stderr,"%g %g\n",x[i],y[i]);
//	}
	fprintf(stderr,"#m=1,S=0\n");
	for(double z=x[0];z<=x[n-1];z+=dz) fprintf(outFileStream_lin,"%g\t %g\t %g\n",z,linterp(n,x,y,z),linterp_integ(n,x,y,z));
//	printf("\n\n");
//	for(double z=x[0];z<=x[n-1];z+=dz) printf(stderr,"%g %g\n",z,linterp_integ(n,x,y,z);
        
fclose(outFileStream_lin);

gsl_interp_free(linear);
return 0;
}

