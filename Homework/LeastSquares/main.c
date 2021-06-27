#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_matrix.h>

#include "extrafuncs.h"
#include "input2Array.h"
#include "gramSchmidt.h"
#include "leastSquares.h"

double funcs(int order, double x){
    switch(order){
        case  0: return 1  ; break;
        case  1: return x  ; break;
        default: return NAN;
    }
}
int main(int argc, char* argv[]){

  int Pts      =  9;
  int Functions    =  2;

  char* InputFile   =  argv[1];
  char* OutputFile  =  argv[2];

  double* xData       =  malloc( Pts*sizeof(double) );
  double* yData       =  malloc( Pts*sizeof(double) );
  double* yDataTrans  =  malloc( Pts*sizeof(double) );
  double* yDev        =  malloc( Pts*sizeof(double) );
  double* yDevTrans   =  malloc( Pts*sizeof(double) );

  input2Array( xData, yData, InputFile );
  compute_devs( Pts, yData, yDev);
  log_featureTransform( Pts, yData, yDataTrans, yDev, yDevTrans );

  gsl_matrix* dataMat    =  gsl_matrix_alloc(Pts, Functions);
  gsl_vector* dataVec    =  gsl_vector_alloc(Pts);
  gsl_vector* coeffsVec  =  gsl_vector_alloc(Functions);
  gsl_matrix* covMat     =  gsl_matrix_alloc(Functions, Functions);


  covMat = leastSquares( Pts, Functions, dataMat, dataVec, coeffsVec, &funcs, xData, yDataTrans, yDevTrans );
  coeffs_exp_featureTransform(coeffsVec, covMat);

  double scale   =  gsl_vector_get(coeffsVec, 0); // a
  double lambda  =  gsl_vector_get(coeffsVec, 1); // -lambda
  write_coeffs(OutputFile, lambda, scale, sqrt(gsl_matrix_get(covMat, 0, 0)), sqrt(gsl_matrix_get(covMat, 1, 1)));

  printf("\n Fit coefficients with errors: \n");
  printf("C_0 (a)      = %lg    +/- %lg   \n", scale,  sqrt(gsl_matrix_get(covMat, 0, 0)));
  printf("C_1 (Lambda) =  %lg +/- %lg \n", lambda, sqrt(gsl_matrix_get(covMat, 1, 1)));

  print_matrix(Functions, covMat, "Cov matrix: ");

  printf("\nhalf-life of ThX     \tt_1/2 = %lg +/- %lg days\n", log(2)/(-lambda), log(2)/(sqrt(gsl_matrix_get(covMat, 1, 1))));
  printf("Should be the same as 224Ra \tt_1/2 = %lg days\n", 3.66                );


  return 0;
}
