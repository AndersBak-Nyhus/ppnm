#ifndef EIGENVALUES_JACOBI_H
#define EIGENVALUES_JACOBI_H

#include <gsl/gsl_matrix.h>

void jacobiMultiply_right(gsl_matrix* matrix, int firstId, int secondId, double angle);
void jacobiMultiply_left( gsl_matrix* matrix, int firstId, int secondId, double angle);
void jacobi_Diag(gsl_matrix* matrix, gsl_matrix* eigVecMat);

#endif
