//
// Created by marc on 3/24/21.
//

#ifndef LEAST_SQUARES_LEASTSQUARES_H
#define LEAST_SQUARES_LEASTSQUARES_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void set_data(int         numOfPts, int         numOfFuncs, double (*fitFuncs)(int, double),
              gsl_matrix* dataMat , gsl_vector* dataVec   ,
              double*     xData   , double*     yData     , double* yDev                   );

gsl_matrix* leastSquares(   int         numOfPts,   int         numOfFuncs,
                            gsl_matrix* dataMat,    gsl_vector* dataVec,
                            gsl_vector* coeffsVec,  double      (*fitFuncs)(int, double),
                            double*     xData    ,  double*     yData                ,  double* yDev  );


#endif //LEAST_SQUARES_LEASTSQUARES_H
