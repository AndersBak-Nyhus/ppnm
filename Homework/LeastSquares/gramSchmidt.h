#ifndef HAVE_GRAMSCHMIDT_DECOMP_H
#define HAVE_GRAMSCHMIDT_DECOMP_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void gramSchmidt_decomp( gsl_matrix* matrixToQR, gsl_matrix* inputTriangularMatrix );
void gramSchmidt_solve(gsl_matrix* orthogonalMatrix, gsl_matrix* triangularMatrix, gsl_vector* rhsVec, gsl_vector* var);
void gramSchmidt_inverse( gsl_matrix* orthogonalMatrix  ,
                          gsl_matrix* triangularMatrix  ,
                          gsl_matrix* inverseMatrix      );
void gramSchmidt_inverseTriang( gsl_matrix* triangularMatrix, gsl_matrix* inverseMatrix );

#endif
