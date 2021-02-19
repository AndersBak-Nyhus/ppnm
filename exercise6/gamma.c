#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_gamma.h>

//double Gamma(double);

int main(){
        double xmin=-2,xmax=2;
        for(double x=xmin;x<=xmax;x+=1.0/8){
                printf("%10g %10g %10g\n",x,tgamma(x),gsl_sf_gamma(x));
        }
return 0;
}
