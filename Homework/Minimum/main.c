#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_vector.h>

#include "minimization.h"

double energy;
double bound_energy;
double maxPt;
double bound_maxPt;

double testFunc (gsl_vector* vals){
    // 1 + (x-a)^2 + (y-b)^2
    double x = gsl_vector_get(vals, 0);
    double y = gsl_vector_get(vals, 1);

    double a = 6;
    double b = 13;

    double funcVal = 1 + (x-a)*(x-a) + (y-b)*(y-b);

    return funcVal;
}

void rosenbrockValley_grad (gsl_vector* vals, gsl_vector* funcVals) {
    double x = gsl_vector_get(vals, 0);
    double y = gsl_vector_get(vals, 1);

    // Set f(x) to be the gradient, f(x) = grad(x)
    gsl_vector_set(funcVals, 0, (-1)*2*(1 - x) + (-2*x)*2*100*(y - x*x));
    gsl_vector_set(funcVals, 1, 2*100*(y - x*x));
}

int main(int argc, char* argv[]){
    // _________________________________________________________________________________________________________________

    printf( "\n################################  " );
    printf( "\n# ------- MINIMIZATION ------- #  " );
    printf( "\n################################\n" );
    printf("\n-- A) Quasi-Newton method with numerical gradient, back-tracking linesearch, and rank-1 update --\n");

    printf("\nTesting the minimization routine, on f(x) = 1 + x² + y² ...\n");
    double a = 6;
    double b = 13;
    double initial_x_val = a - 2;
    double initial_y_val = b - 3;

    int dims = 2;
    double tolerance = 1e-5;
    gsl_vector* minimum = gsl_vector_alloc(dims);
    gsl_vector_set(minimum, 0, initial_x_val );
    gsl_vector_set(minimum, 0, initial_y_val );

    printf("The initial value is at (x, y) = (%g, %g)\n", initial_x_val, initial_y_val);
    quasi_newtonMethod(testFunc, minimum, tolerance);
    printf("The found  minima is at (x, y) = (%g, %g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));
    printf("The actual minima is at (x, y) = (%g, %g)\n", a, b);

    return 0;
}
