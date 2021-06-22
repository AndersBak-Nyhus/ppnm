#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "komplex.h"

void komplex_print (char *s, komplex a) {
	printf ("%s (%g,%g)\n", s, a.re, a.im);
}

void komplex_set (komplex* z, double x, double y) {
	(*z).re = x;
	(*z).im = y;
}

komplex komplex_new (double x, double y) {
	komplex z = { x, y };
	return z;
}

komplex komplex_add (komplex a, komplex b) {
	komplex result = { a.re + b.re , a.im + b.im };
	return result;
}

komplex komplex_sub (komplex a, komplex b) {
	komplex result = { a.re - b.re , a.im - b.im };
	return result;
}

int komplex_equal (komplex a, komplex b, double acc, double eps){
	 if (
		 	 (    (fabs( a.re - b.re ) < acc)
	      &&  (fabs( a.im - b.im ) < acc)
			 )
	 	    ||
			 (   ((fabs( a.re - b.re ) / ( fabs(a.re) + fabs(b.re) )) < (eps / 2))
	      && ((fabs( a.im - b.im ) / ( fabs(a.re) + fabs(b.re) )) < (eps / 2))
			 )
		  )  { return 1; }
	  else { return 0; }
}

komplex komplex_mul (komplex a, komplex b) {
	komplex result = komplex_new( a.re * b.re - a.im * b.im, a.im * b.re + a.re * b.im );
	return result;
}

komplex komplex_div (komplex a, komplex b){
	komplex result = komplex_new( (a.re*b.re + a.im*b.im)/(b.re*b.re + b.im*b.im), (a.im*b.re - a.re*b.im)/(b.re*b.re + b.im*b.im) );
	return result;
}

komplex komplex_conjugate(komplex z){
	z.im = -(z.im);
	return z;
}

double  komplex_abs      (komplex z){
	double abs = sqrt(z.re*z.re + z.im*z.im);
	return abs;
}

komplex komplex_exp      (komplex z){
	// The exponential of a complex number will be
	// exp( z.re + i*z.im ) = exp(z.re)*exp(i*z.im)
	// where exp(i*z.im) = cos(z.im)+i*sin(z.im)
	// hence:
	komplex expRe   =  komplex_new( exp( z.re ), 0);
	komplex expIm  =  komplex_new( cos(z.im), sin(z.im) );
	return komplex_mul(expRe, expIm);
}

komplex komplex_sin      (komplex z){
	komplex m_I     =  komplex_new(0, 1);
	komplex expArg  =  komplex_mul(m_I, z);
	komplex mexpArg =  komplex_mul(komplex_new(-1,0), expArg);
	komplex num 		=  komplex_sub(komplex_exp(expArg), komplex_exp(mexpArg));
	komplex den 		=  komplex_mul(komplex_new(2,0), m_I);
	return komplex_div(num, den);
}
komplex komplex_cos      (komplex z){
	komplex m_I 		=	 komplex_new(0, 1);
	komplex expArg 	=	 komplex_mul(m_I, z);
	komplex mexpArg =  komplex_mul(komplex_new(-1,0), expArg);
	komplex num 		=	 komplex_add(komplex_exp(expArg), komplex_exp(mexpArg));
	komplex den 		=	 komplex_new(2,0);
	return komplex_div(num, den);
}

komplex komplex_sqrt     (komplex z){
	return komplex_new(sqrt( (komplex_abs(z)+z.re)/2), (z.im/abs(z.im))*sqrt((komplex_abs(z)-z.re)/2) );
}
