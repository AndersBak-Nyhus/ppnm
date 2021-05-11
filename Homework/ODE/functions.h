//
// Created by marc on 4/12/21.
//

#ifndef CLION_TESTING_FUNCTIONS_H
#define CLION_TESTING_FUNCTIONS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void harmonicFunc   (double var, gsl_vector* funcVal, gsl_vector* funcDeriv);
void SIRmodel       (double var, gsl_vector* funcVal, gsl_vector* funcDeriv);
void SIRmodel2      (double var, gsl_vector* funcVal, gsl_vector* funcDeriv);
void threeBodyProb  (double var, gsl_vector* funcVal, gsl_vector* funcDeriv);

#endif //CLION_TESTING_FUNCTIONS_H
