#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "gramSchmidt.h"
#include "extraFuncs.h"
#include "timer.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

void runtime(int Rep, int startRep, char* outputFilename, char* gsl_outputFilename, unsigned int* seed);

int main (int argc, char* argv[]){

  unsigned int seed   =   time(NULL) ;

  int Rep = 2.5e2;
  int startRep = 200;

  int Rows = 3;
	int Cols = 2;
//allocate memory
  gsl_matrix* ortgnlMatTall    =  gsl_matrix_alloc(Rows, Cols);
  gsl_matrix* ortgnlMatSquare  =  gsl_matrix_alloc(Rows, Rows);
  gsl_matrix* testMatTall      =  gsl_matrix_alloc(Rows, Cols);
  gsl_matrix* testMatSquare    =  gsl_matrix_alloc(Rows, Rows);
  gsl_matrix* triangMatTall    =  gsl_matrix_alloc(Cols, Cols);
  gsl_matrix* triangMatSquare  =  gsl_matrix_alloc(Rows, Rows);
  gsl_vector* RHSvecTall       =  gsl_vector_alloc(Rows);
  gsl_vector* RHSvecSquare     =  gsl_vector_alloc(Rows);
  gsl_vector* varVec           =  gsl_vector_alloc(Rows);
  gsl_vector* varVecTall       =  gsl_vector_alloc(Cols);
  gsl_matrix* checkOrtgnl      =  gsl_matrix_alloc(Cols, Cols);
  gsl_matrix* checkRes         =  gsl_matrix_alloc(Rows, Cols);
  gsl_vector* checkRHS         =  gsl_vector_alloc(Rows);
  gsl_matrix* inverseMatrix    =  gsl_matrix_alloc(Rows, Rows);
  gsl_matrix* checkInverse     =  gsl_matrix_alloc(Rows, Rows);

//Set data in matrix
	set_data(testMatTall, RHSvecTall, &seed);
  set_data(testMatSquare, RHSvecSquare, &seed);
  gsl_matrix_memcpy(ortgnlMatTall, testMatTall);

//decomposition
  gramSchmidt_decomp(ortgnlMatTall, triangMatTall);


//print matrix
  print_matrix(Cols, triangMatTall, "Part A)\nCheck of upper triangular matrix:");

//check Q^T*Q
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, ortgnlMatTall, ortgnlMatTall, 0, checkOrtgnl);
  print_matrix(Cols, checkOrtgnl, "check that Q^T * Q = 1:");

//check QR
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, ortgnlMatTall, triangMatTall, 0, checkRes);
  print_matrix(Rows, checkRes, "QR:");
  print_matrix(Rows, testMatTall, "Tall matrix:");

//solver
  print_matrix(Rows, testMatSquare, "Square matrix for solver:");
  gsl_matrix_memcpy(ortgnlMatSquare, testMatSquare);
  gramSchmidt_decomp(ortgnlMatSquare, triangMatSquare);
  gramSchmidt_solve(ortgnlMatSquare, triangMatSquare, RHSvecSquare, varVec);
  gsl_blas_dgemv(CblasNoTrans, 1, testMatSquare, varVec, 0, checkRHS);
  vector_print( "\nFrom gramSchmidt_solve             Ax = ", checkRHS);
  vector_print( "\nShould be equal to right-hand side  b = ", RHSvecSquare);


//Inverse
  gramSchmidt_inverse(ortgnlMatSquare, triangMatSquare, inverseMatrix);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, testMatSquare, inverseMatrix, 0, checkInverse);
  print_matrix(Rows, checkInverse, "\nPart B)\nChecking A*A^{-1} = I :");

//free allocated memory
  gsl_matrix_free (ortgnlMatTall);
  gsl_matrix_free (ortgnlMatSquare);
	gsl_matrix_free (testMatTall);
	gsl_matrix_free (testMatSquare);
	gsl_matrix_free (triangMatTall);
	gsl_matrix_free (triangMatSquare);
	gsl_vector_free (RHSvecTall);
	gsl_vector_free (RHSvecSquare);
	gsl_vector_free (varVec);
  gsl_matrix_free (checkOrtgnl);
  gsl_matrix_free (checkRes);
	gsl_vector_free (checkRHS);
  gsl_matrix_free (inverseMatrix);
  gsl_matrix_free (checkInverse);

//check runtime of QR
  runtime(Rep, startRep, argv[1], argv[2], &seed);

  return 0;
}
