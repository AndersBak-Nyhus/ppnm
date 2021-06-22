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

  int NumberOfPts      =  9;
  int NumberOfFuncs    =  2;

  char* InputFile   =  argv[1];
  char* OutputFile  =  argv[2];

  double* xData       =  malloc( NumberOfPts*sizeof(double) );
  double* yData       =  malloc( NumberOfPts*sizeof(double) );
  double* yDataTrans  =  malloc( NumberOfPts*sizeof(double) );
  double* yDev        =  malloc( NumberOfPts*sizeof(double) );
  double* yDevTrans   =  malloc( NumberOfPts*sizeof(double) );

  input2Array( xData, yData, InputFile );
  compute_deviations( NumberOfPts, yData, yDev);
  log_featureTransform( NumberOfPts, yData, yDataTrans, yDev, yDevTrans );

  gsl_matrix* dataMat    =  gsl_matrix_alloc(NumberOfPts,   NumberOfFuncs);
  gsl_vector* dataVec    =  gsl_vector_alloc(NumberOfPts              );
  gsl_vector* coeffsVec  =  gsl_vector_alloc(NumberOfFuncs            );
  gsl_matrix* covMat     =  gsl_matrix_alloc(NumberOfFuncs, NumberOfFuncs);


  covMat = leastSquares( NumberOfPts, NumberOfFuncs, dataMat, dataVec, coeffsVec, &funcs, xData, yDataTrans, yDevTrans );
  coeffs_exp_featureTransform(coeffsVec, covMat);

  double scale   =  gsl_vector_get(coeffsVec, 0); // a
  double lambda  =  gsl_vector_get(coeffsVec, 1); // -lambda
  write_coeffs(OutputFile, lambda, scale, sqrt(gsl_matrix_get(covMat, 0, 0)), sqrt(gsl_matrix_get(covMat, 1, 1)));

  printf("\n Fit coefficients: \n");
  printf("C_0 (a)      = %lg    +/- %lg   \n", scale,  sqrt(gsl_matrix_get(covMat, 0, 0)));
  printf("C_1 (Lambda) =  %lg +/- %lg \n", lambda, sqrt(gsl_matrix_get(covMat, 1, 1)));

  print_matrix(NumberOfFuncs, covMat, "Covariance matrix: ");

  printf("\nhalf-life of ThX     \tt_1/2 = %lg +/- %lg days\n", log(2)/(-lambda), log(2)/(sqrt(gsl_matrix_get(covMat, 1, 1))));
  printf("Should be the same as 224Ra \tt_1/2 = %lg days\n", 3.66                );


  return 0;
}
