#ifndef HAVE_UTILITIES_H
#define HAVE_UTILITIES_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

double randomNumber    ( unsigned int *seed );
void   vector_print    (char* string, gsl_vector* vector);
void   set_data_tall   (gsl_matrix* testMatTall, gsl_vector* RHSvecTall, unsigned int *LHSseed, unsigned int *RHSseed);
void   set_data_square (gsl_matrix* testMatSquare, gsl_vector* RHSvecSquare, unsigned int *LHSseed, unsigned int *RHSseed);
void   print_matrix    (int numOfRows, gsl_matrix* matrixToPrint, char* string );
void   compute_devs(int numOfPts, double* yData, double* yDev);
void   log_featureTransform( int numOfPts, double* yData, double* yDataTrans, double* yDev, double* yDevTrans );
void   coeffs_exp_featureTransform(gsl_vector* coeffsVec, gsl_matrix* covMat);
void   write_coeffs(char* outputFilename, double lambda, double scale, double scale_err, double lambda_err);

#endif
