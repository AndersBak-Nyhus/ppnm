#ifndef HAVE_MINIMIZATION_H
#define HAVE_MINIMIZATION_H

#include <gsl/gsl_vector.h>

void numeric_gradient (double func(gsl_vector*), gsl_vector* minimum, gsl_vector* gradient);
void quasi_newtonMethod( double func(gsl_vector*), gsl_vector* minimum, double tolerance);

#endif
