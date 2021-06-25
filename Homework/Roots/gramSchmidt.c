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

  //Modified gram-schmidt orthogonalization
  for (int colId = 0; colId < Cols; colId++ ){

    gsl_vector* col   =   gsl_vector_alloc(Rows);
    *col              =   (gsl_matrix_column( matrixToQR, colId )).vector;
    double colNorm    =   gsl_blas_dnrm2(col);                                  

    gsl_matrix_set( inputTriangularMatrix, colId, colId, colNorm );             


   
    gsl_vector* orthgnMatCol  =   gsl_vector_alloc(Rows);                  

    gsl_vector_memcpy(orthgnMatCol, col);                                       
    
    gsl_vector_scale(orthgnMatCol, 1.0/colNorm);                                
    
    gsl_matrix_set_col(matrixToQR, colId, orthgnMatCol);                        

   
    for (int nextColId = colId + 1; nextColId < Cols; nextColId++ ){           
    

      gsl_vector* nextCol   =   gsl_vector_alloc(Rows);
      *nextCol              =   (gsl_matrix_column( matrixToQR, nextColId )).vector;

      double triangMatElement;                                                      
      
      gsl_blas_ddot(orthgnMatCol, nextCol, &triangMatElement);                      
      
      gsl_matrix_set( inputTriangularMatrix, colId, nextColId, triangMatElement );  
      
      gsl_vector* ortColScaled = gsl_vector_alloc(Rows);                       
      
      gsl_vector_memcpy(ortColScaled, orthgnMatCol);                                
      
      gsl_vector_scale(ortColScaled, triangMatElement);                             
      
      gsl_vector_sub(nextCol, ortColScaled);                                        
                                                                                    

      gsl_matrix_set_col(matrixToQR, nextColId, nextCol);                           
                                                                                    
                                                                                    
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
