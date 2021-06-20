#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "gramSchmidt.h"
#include "extraFuncs.h"
#include "diffClock.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

void test_runtime(int numOfReps, int startRep, char* my_outputFilename, char* gsl_outputFilename, unsigned int* seed);

int main (int argc, char* argv[]){
  if (argc < 2){
    fprintf(stdout, "No input arguments taken.");
    exit(-1);
  }
  unsigned int seed   =   time(NULL) ;

  int numOfReps = 2.5e2;
  int startRep = 200;

  int numOfRows = 3;
	int numOfCols = 2;

  gsl_matrix* ortgnlMatTall    =  gsl_matrix_alloc(numOfRows, numOfCols);
  gsl_matrix* ortgnlMatSquare  =  gsl_matrix_alloc(numOfRows, numOfRows);
  gsl_matrix* testMatTall      =  gsl_matrix_alloc(numOfRows, numOfCols);
  gsl_matrix* testMatSquare    =  gsl_matrix_alloc(numOfRows, numOfRows);
  gsl_matrix* triangMatTall    =  gsl_matrix_alloc(numOfCols, numOfCols);
  gsl_matrix* triangMatSquare  =  gsl_matrix_alloc(numOfRows, numOfRows);
  gsl_vector* RHSvecTall       =  gsl_vector_alloc(numOfRows           );
  gsl_vector* RHSvecSquare     =  gsl_vector_alloc(numOfRows           );
  gsl_vector* varVec           =  gsl_vector_alloc(numOfRows           );
  gsl_vector* varVecTall       =  gsl_vector_alloc(numOfCols           );
  gsl_matrix* checkOrtgnl      =  gsl_matrix_alloc(numOfCols, numOfCols);
  gsl_matrix* checkRes         =  gsl_matrix_alloc(numOfRows, numOfCols);
  gsl_vector* checkRHS         =  gsl_vector_alloc(numOfRows           );
  gsl_matrix* inverseMatrix    =  gsl_matrix_alloc(numOfRows, numOfRows);
  gsl_matrix* checkInverse     =  gsl_matrix_alloc(numOfRows, numOfRows);

	set_data(testMatTall, RHSvecTall, &seed);
  set_data(testMatSquare, RHSvecSquare, &seed);
  gsl_matrix_memcpy(ortgnlMatTall, testMatTall);

  gramSchmidt_decomp(ortgnlMatTall, triangMatTall);

  print_matrix(numOfCols, triangMatTall, "Triangular matrix:");

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, ortgnlMatTall, ortgnlMatTall, 0, checkOrtgnl);
  print_matrix(numOfCols, checkOrtgnl, "Q^T * Q:");

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, ortgnlMatTall, triangMatTall, 0, checkRes);
  print_matrix(numOfRows, checkRes, "QR :");
  print_matrix(numOfRows, testMatTall, "Tall test matrix:");

  print_matrix(numOfRows, testMatSquare, "Square test matrix for solver:");
  gsl_matrix_memcpy(ortgnlMatSquare, testMatSquare);
  gramSchmidt_decomp(ortgnlMatSquare, triangMatSquare);
  gramSchmidt_solve(ortgnlMatSquare, triangMatSquare, RHSvecSquare, varVec);
  gsl_blas_dgemv(CblasNoTrans, 1, testMatSquare, varVec, 0, checkRHS);
  vector_print( "\nFrom gramSchmidt_solve             Ax = ", checkRHS);
  vector_print( "\nShould be equal to right-hand side  b = ", RHSvecSquare);

  gramSchmidt_inverse(ortgnlMatSquare, triangMatSquare, inverseMatrix);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, testMatSquare, inverseMatrix, 0, checkInverse);
  print_matrix(numOfRows, checkInverse, "\nChecking A*A^{-1} = I :");

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

  test_runtime(numOfReps, startRep, argv[1], argv[2], &seed);

  return 0;
}
