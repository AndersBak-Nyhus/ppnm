#ifndef CLION_TESTING_RK_ODE_H
#define CLION_TESTING_RK_ODE_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void rkstep12 ( void        (*func)(double, gsl_vector*, gsl_vector*) ,
                double      var                                       ,
                gsl_vector* funcVal                                   ,
                double      step                                      ,
                gsl_vector* funcStep                                  ,
                gsl_vector* err                                         );
void rkstep45(  void        (*func)(double, gsl_vector*, gsl_vector*) ,
                double      var                                       ,
                gsl_vector* funcVal                                   ,
                double      step                                      ,
                gsl_vector* funcStep                                  ,
                gsl_vector* err                                         );
void rkdriver(  void        (*func)(double, gsl_vector*, gsl_vector*) ,
                double      leftEndpt                                 ,
                gsl_vector* funcValLeft                               ,
                double      rightEndpt                                ,
                gsl_vector* funcValRight                              ,
                double      step                                      ,
                double      absAcc                                    ,
                double      relAcc                                    ,
                FILE*       path2File                                   );

#endif 
