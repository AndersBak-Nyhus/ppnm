#include<math.h>
#include<stdio.h>
double ex(double);
int main(){
	for(double x=1./16;x<=8;x+=1./16)
		printf("%g %g %g\n",x,ex(x),exp(x));

return 0;
}
