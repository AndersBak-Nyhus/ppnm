#ifndef HAVE_UTILITIES_H
#define HAVE_UTILITIES_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

double randomNumber    ( unsigned int *seed );
void   vector_print    (char* string, gsl_vector* vector);
void   set_data        (gsl_matrix* tallMatrix, gsl_vector* tallRHSVector, unsigned int *seed);
void   print_matrix    (int numOfRows, gsl_matrix* matrixToPrint, char* string );

#endif
