#include <assert.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "gramSchmidt.h"
#include "backsub.h"

void gramSchmidt_decomp( gsl_matrix* matrixToQR, gsl_matrix* inputTriangularMatrix ){


  int Rows   =  (int) matrixToQR -> size1;
  int Cols   =  (int) matrixToQR -> size2;

  assert( Rows >= Cols );

  // Do modified gram-schmidt orthogonalization
  for (int Colnum = 0; Colnum < Cols; Colnum++ ){

    // Begin by normalizing the i'th (Colnum) column of A (matrixToQR), a_i.
    gsl_vector* col   =   gsl_vector_alloc(Rows);
    *col              =   (gsl_matrix_column( matrixToQR, Colnum )).vector;
    double colNorm    =   gsl_blas_dnrm2(col);                                  //  Compute euclidean R^2-norm

    gsl_matrix_set( inputTriangularMatrix, Colnum, Colnum, colNorm );             //  Set diagonal element of triangular matrix equal to this norm


    // Compute the i'th (Colnum) column of the orthogonal matrix Q
    gsl_vector* orthgnMatCol  =   gsl_vector_alloc(Rows);                  //  gsl_vector to hold column

    gsl_vector_memcpy(orthgnMatCol, col);                                       //  Copy contents of col into orthgnMatCol
    gsl_vector_scale(orthgnMatCol, 1.0/colNorm);                                //  Scale orthognMatCol, dividing each element by the norm
    gsl_matrix_set_col(matrixToQR, Colnum, orthgnMatCol);                        //  Set the i'th (Colnum) column of A (matrixToQR) equal to the orthogonal column

    // Make the remaining columns of the matrix orthogonal to the new normalized column
    for (int nextColnum = Colnum + 1; nextColnum < Cols; nextColnum++ ){           //  Loop over remaining columns j (nextColnum)

      gsl_vector* nextCol   =   gsl_vector_alloc(Rows);
      *nextCol              =   (gsl_matrix_column( matrixToQR, nextColnum )).vector;

      double triangMatElement;                                                      //  Matrix element of upper triangular matrix, R_ij
      gsl_blas_ddot(orthgnMatCol, nextCol, &triangMatElement);                      //  Compute the matrix element R_ij = q_i^T * a_j
      gsl_matrix_set( inputTriangularMatrix, Colnum, nextColnum, triangMatElement );  //  Set the ij'th element of the upper triangular matrix to be R_ij (triangMatElement)
      gsl_vector* ortColScaled = gsl_vector_alloc(Rows);                       //  Vector to hold the scaled i'th (Colnum) orthognal vector
      gsl_vector_memcpy(ortColScaled, orthgnMatCol);                                //  Copy q_i, the i'th (Colnum) orthogonal column into this new vector
      gsl_vector_scale(ortColScaled, triangMatElement);                             //  Compute q_i*R_ij
      gsl_vector_sub(nextCol, ortColScaled);                                        //  Compute orthogonal complement by subtracting the component along q_i, q_i*R_ij, from a_j,
                                                                                    //  that is, do exactly a_j = a_j - q_i*R_ij

      gsl_matrix_set_col(matrixToQR, nextColnum, nextCol);                           //  Go back and set the j'th (nextColnum) column of A (matrixToQR)
                                                                                    //  equal to the column from before. This makes all the reamining columns of
                                                                                    //  A (matrixToQR), from j = i + 1,, orthogonal to the i'th orhogonal column q_i
      gsl_vector_free(nextCol);
      gsl_vector_free(ortColScaled);
    }
    gsl_vector_free(col);
    gsl_vector_free(orthgnMatCol);
  }
}


void gramSchmidt_solve( gsl_matrix* orthogonalMatrix  ,
                        gsl_matrix* triangularMatrix  ,
                        gsl_vector* rhsVec            ,
                        gsl_vector* var                 ){


  gsl_blas_dgemv(CblasTrans, 1.0, orthogonalMatrix, rhsVec, 0.0, var);
  backsub(triangularMatrix, var);

}


void gramSchmidt_inverse( gsl_matrix* orthogonalMatrix  ,
                          gsl_matrix* triangularMatrix  ,
                          gsl_matrix* inverseMatrix      ){


  int Dims = (triangularMatrix -> size1);

  gsl_matrix* triangInverse  =  gsl_matrix_alloc(Dims, Dims);
  gsl_vector* unitVec        =  gsl_vector_alloc(Dims);

  for (int dim = 0; dim < Dims; dim++ ){
    gsl_vector_set_basis(unitVec, dim);
    backsub(triangularMatrix, unitVec);
    gsl_matrix_set_col(triangInverse, dim, unitVec);
  }

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, triangInverse, orthogonalMatrix, 0, inverseMatrix);

  gsl_matrix_free(triangInverse);
  gsl_vector_free(unitVec);
}
