#ifndef LEAST_SQUARES_LEASTSQUARES_H
#define LEAST_SQUARES_LEASTSQUARES_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void set_data(int Pts, int Functions, double (*fitFuncs)(int, double),
              gsl_matrix* dataMat , gsl_vector* dataVec   ,
              double* xData, double* yData, double* yDev);

gsl_matrix* leastSquares(int Pts,   int Functions,
                            gsl_matrix* dataMat,    gsl_vector* dataVec,
                            gsl_vector* coeffsVec,  double (*fitFuncs)(int, double),
                            double* xData, double* yData , double* yDev);


#endif
