//
// Created by marc on 4/12/21.
//

#include "functions.h"
#include <math.h>
#include <gsl/gsl_vector.h>

void harmonicFunc(double var, gsl_vector* funcVal, gsl_vector* funcDeriv){
    /*
     *  A testing function, the harmonic oscillator
     *
     *      d^2 u/dt = -u
     *
     *  written as a first order planar equation
     *
     *      dy/dt = f(t,y)
     *
     *  which simply takes a variable (var) representing t, and a vector y,
     *  and then fills out funcDeriv, a vector, that represents the RHS
     *  of the above ODE. For use with the RK-ODE.h methods.
     *
     */

    gsl_vector_set(funcDeriv, 0,   gsl_vector_get(funcVal, 1));
    gsl_vector_set(funcDeriv, 1, - gsl_vector_get(funcVal, 0));
}

void SIRmodel(double var, gsl_vector* funcVal, gsl_vector* funcDeriv){
    /*
     *  SIR-model function implementation as a six dimensional
     *  planar equation to use with the ODE solver from RK-ODE.h
     *  which implements a runge-kutta-method ode-solver.
     *
     */

    double populationSize  =  5808180;  // Populations size of Denmark
    double contactTime     =  2.5;      // Contact time in days
    double recoveryTime    =  25;       // Recovery time in days

    double susceptible  =  gsl_vector_get(funcVal, 0);
    double infected     =  gsl_vector_get(funcVal, 1);

    double diffSusceptible_RHS  = - (infected * susceptible) / (populationSize * contactTime);
    double diffRemoved_RHS      =   infected / recoveryTime;
    double diffInfected_RHS     = - diffSusceptible_RHS - diffRemoved_RHS;

    gsl_vector_set(funcDeriv, 0, diffSusceptible_RHS);
    gsl_vector_set(funcDeriv, 1, diffInfected_RHS);
    gsl_vector_set(funcDeriv, 2, diffRemoved_RHS);
}

void SIRmodel2(double var, gsl_vector* funcVal, gsl_vector* funcDeriv){
    /*
     *  This function is the same model as the one just above
     *  instead here T_c has been double to study the effect
     *  of increasing values of this parameter.
     *
     */

    double populationSize  =  5808180;  // Populations size of Denmark
    double contactTime     =  5;        // Contact time, HERE DOUBLED
    double recoveryTime    =  25;       // Recovery time in days

    double susceptible  =  gsl_vector_get(funcVal, 0);
    double infected     =  gsl_vector_get(funcVal, 1);

    double diffSusceptible_RHS  = - (infected * susceptible) / (populationSize * contactTime);
    double diffRemoved_RHS      =   infected / recoveryTime;
    double diffInfected_RHS     = - diffSusceptible_RHS - diffRemoved_RHS;

    gsl_vector_set(funcDeriv, 0, diffSusceptible_RHS);
    gsl_vector_set(funcDeriv, 1, diffInfected_RHS);
    gsl_vector_set(funcDeriv, 2, diffRemoved_RHS);
}

void threeBodyProb  (double var, gsl_vector* funcVal, gsl_vector* funcDeriv){
    /*
     * The three body problem function for use with the ODE-solver methods
     * found in RK-ODE.h, here implemented for the special case of equal
     * masses of the three particles.
     *
     */
    double gravitationalConst   =   1;
    double mass_1               =   1;
    double mass_2               =   mass_1;
    double mass_3               =   mass_1;

    double pos_x_1  =  gsl_vector_get(funcVal, 0);
    double pos_y_1  =  gsl_vector_get(funcVal, 1);
    double pos_x_2  =  gsl_vector_get(funcVal, 2);
    double pos_y_2  =  gsl_vector_get(funcVal, 3);
    double pos_x_3  =  gsl_vector_get(funcVal, 4);
    double pos_y_3  =  gsl_vector_get(funcVal, 5);

    gsl_vector_set(funcDeriv, 0, gsl_vector_get(funcVal, 6 ));
    gsl_vector_set(funcDeriv, 1, gsl_vector_get(funcVal, 7 ));
    gsl_vector_set(funcDeriv, 2, gsl_vector_get(funcVal, 8 ));
    gsl_vector_set(funcDeriv, 3, gsl_vector_get(funcVal, 9 ));
    gsl_vector_set(funcDeriv, 4, gsl_vector_get(funcVal, 10));
    gsl_vector_set(funcDeriv, 5, gsl_vector_get(funcVal, 11));


    double dist_12 = sqrt(pow((pos_x_1 - pos_x_2), 2) + pow((pos_y_1 - pos_y_2), 2));
    double dist_23 = sqrt(pow((pos_x_2 - pos_x_3), 2) + pow((pos_y_2 - pos_y_3), 2));
    double dist_13 = sqrt(pow((pos_x_3 - pos_x_1), 2) + pow((pos_y_3 - pos_y_1), 2));

    gsl_vector_set(funcDeriv, 6,  (pos_x_2 - pos_x_1) / pow(dist_12, 3) + (pos_x_3 - pos_x_1) / pow(dist_13, 3));
    gsl_vector_set(funcDeriv, 7,  (pos_y_2 - pos_y_1) / pow(dist_12, 3) + (pos_y_3 - pos_y_1) / pow(dist_13, 3));
    gsl_vector_set(funcDeriv, 8,  (pos_x_3 - pos_x_2) / pow(dist_23, 3) + (pos_x_1 - pos_x_2) / pow(dist_12, 3));
    gsl_vector_set(funcDeriv, 9,  (pos_y_3 - pos_y_2) / pow(dist_23, 3) + (pos_y_1 - pos_y_2) / pow(dist_12, 3));
    gsl_vector_set(funcDeriv, 10, (pos_x_2 - pos_x_3) / pow(dist_23, 3) + (pos_x_1 - pos_x_3) / pow(dist_13, 3));
    gsl_vector_set(funcDeriv, 11, (pos_y_2 - pos_y_3) / pow(dist_23, 3) + (pos_y_1 - pos_y_3) / pow(dist_13, 3));
}