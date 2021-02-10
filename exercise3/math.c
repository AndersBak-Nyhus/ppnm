#include<stdio.h>
#include<math.h>
#include<complex.h>

int main(){
	//calculating gamma(5)
	double g=tgamma(5);
	printf("gamma(5)=%g\n",g);

	//calculating J_1(0.5)
	double b=j1(0.5);
	printf("J1(0.5)=%g\n",b);

	//calculating sqrt(-2)
	complex s;
	s=csqrt(-2);
	printf("sqrt(-2)=%g+I%g\n",creal(s),cimag(s));

	//calculating exp(i*pi)
	complex ss;
	ss=cexp(I*M_PI);
	printf("exp(i*pi)=%g+I%g\n",creal(ss),cimag(ss));

	//calculating exp(i)
	complex z;
	z=cexp(I);
	printf("exp(i)=%g+I%g\n",creal(z),cimag(z));

	//calculating i^(e)
	complex zz;
	zz=cpow(I,M_E);
	printf("i^e=%g+%gI\n",creal(zz),cimag(zz));

	//calculating i^i
	complex r;
	r=cpow(I,I);
	printf("i^i=%g+%gI\n",creal(r),cimag(r));

	//significant digits
	float x_float=1.f/9;
	printf("float\n");
	printf("1/9=%.25g\n",x_float);

	double x_double=1./9;
	printf("double\n");
	printf("1/9=%.25lg\n",x_double);

	long double x_long_double=1.L/9;
	printf("long double\n");
	printf("1/9=%.25Lg\n",x_long_double);
return 0;
}
