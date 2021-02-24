#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

double f(double x,void* params){
	//double z =*(double*)params;
	double f = log(x)/sqrt(x);
return f;
}

double exA(){
	gsl_function F;
	F.function = &f;
	//F.params = (void*)&z;
	int limit = 999;
	gsl_integration_workspace* w;
	w = gsl_integration_workspace_alloc (limit);
	double a=0,b=1,acc=1e-6,eps=1e-6,result,error;
	gsl_integration_qags(&F,a,b,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}

int main(){
	printf("nummerical integration of ln(x)/sqrt(x) from 0 to 1=%10g\n",exA());
return 0;
}
