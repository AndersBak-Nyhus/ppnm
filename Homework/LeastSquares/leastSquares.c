#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "leastSquares.h"
#include "gramSchmidt.h"
#include "extrafuncs.h"

void set_data(int  Pts, int Functions, double (*fitFuncs)(int, double),
              gsl_matrix* dataMat , gsl_vector* dataVec,
              double* x, double* y, double* yDev){
    for (int row = 0; row < Pts; row++){
        for (int col = 0; col < Functions; col++){
            gsl_matrix_set(dataMat, row, col, (fitFuncs(col, x[row]))/yDev[row]);
        }
        gsl_vector_set(dataVec, row, y[row]/yDev[row]);
    }
}

gsl_matrix* leastSquares(  int Pts, int Functions,
                           gsl_matrix* dataMat,    gsl_vector* dataVec,
                           gsl_vector* coeffsVec,  double      (*fitFuncs)(int, double),
                           double* x, double* y,  double* yDev  ){

    set_data(Pts, Functions, fitFuncs, dataMat, dataVec, x, y, yDev);

    gsl_matrix* ortgMat  =  gsl_matrix_alloc(Pts, Functions); 
    gsl_matrix* triangMat =  gsl_matrix_alloc(Functions, Functions); 
    
    gsl_matrix* triangMatInverse  =  gsl_matrix_alloc(Functions, Functions);
    gsl_vector* tmpRes =  gsl_vector_alloc(Functions);

    gsl_matrix_memcpy(ortgMat, dataMat);
    gramSchmidt_decomp(ortgMat, triangMat);
    gramSchmidt_solve(ortgMat, triangMat, dataVec, coeffsVec);

    
    gsl_matrix* covMat  =  gsl_matrix_alloc(Functions, Functions);
    gramSchmidt_inverseTriang( triangMat, triangMatInverse );
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, triangMatInverse, triangMatInverse,  0, covMat);

    return covMat;
}
