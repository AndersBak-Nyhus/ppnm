#ifndef MONTE_CARLO_MONTECARLO_H
#define MONTE_CARLO_MONTECARLO_H

double van_der_corput_sequence(int id , int base);
void halton_sequence(int id, int dim, double *vec);
void halton_sequence_secondary(int id, int dim, double *vec);
void lattice(int d, double* x);
void rand_num (int dim , const double *lowerBound , const double *upperBound , double *nums );
void rand_num_halton_corput (int id, int dim , const double *lowerBound , const double *upperBound , double *vec );
void rand_num_halton_corput_secondary (int id, int dim , const double *lowerBound , const double *upperBound , double *vec );
void plain_montecarlo (int dim, double *lowerBound , double *upperBound, double func(double* nums) , int numOfPts, double* result , double* error );
void plain_montecarlo_quasi (int dim, double *lowerBound , double *upperBound, double func(double* nums) , int numOfPts, double* result , double* error );
double plain_montecarlo_stratifiedSampling(  int     dim                         ,
                                             double  func(double* vec)           ,
                                             double* lowerBound                  ,
                                             double* upperBound                  ,
                                             double  absAcc                      ,
                                             double  relAcc                      ,
                                             int     numOfRecalls                ,
                                             double  meanRecall                      );
#endif //MONTE_CARLO_MONTECARLO_H
