//
// Created by marc on 3/26/21.
//

#ifndef EIGENVALUES_JACOBI_H
#define EIGENVALUES_JACOBI_H

#include <gsl/gsl_matrix.h>

void jacobiMultiply_right(gsl_matrix* matrix, int firstId, int secondId, double angle);
void jacobiMultiply_left( gsl_matrix* matrix, int firstId, int secondId, double angle);
void jacobiDiag(gsl_matrix* matrix, gsl_matrix* eigVecMat);

#endif //EIGENVALUES_JACOBI_H
