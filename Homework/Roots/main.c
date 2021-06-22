#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "gramSchmidt.h"
#include "rootfinding.h"
#include "RK-ODE.h"

double energy;
double bound_energy;
double maxPt;
double bound_maxPt;

void testFunc (gsl_vector* vals, gsl_vector* funcVals){
    // 1 + (x-a)^2 + (y-b)^2
    double x = gsl_vector_get(vals, 0);
    double y = gsl_vector_get(vals, 1);

    double a = 6;
    double b = 13;

    // Set f(x) to be the gradient, f(x) = grad(x)
    gsl_vector_set(funcVals, 0, 2*x*(x-a));
    gsl_vector_set(funcVals, 1, 2*y*(y-b));

}

void rosenbrockValley_grad (gsl_vector* vals, gsl_vector* funcVals) {
    double x = gsl_vector_get(vals, 0);
    double y = gsl_vector_get(vals, 1);

    // Set f(x) to be the gradient, f(x) = grad(x)
    gsl_vector_set(funcVals, 0, (-1)*2*(1 - x) + (-2*x)*2*100*(y - x*x));
    gsl_vector_set(funcVals, 1, 2*100*(y - x*x));
}

void schrodingerEq (double var, gsl_vector* funcVal, gsl_vector* funcDeriv){
    double thisFuncVal = gsl_vector_get(funcVal, 0);
    double firstDeriv = gsl_vector_get(funcVal, 1);
    double secondDeriv = (-2)*(1.0/var + energy)*thisFuncVal;

    gsl_vector_set(funcDeriv, 0, firstDeriv);
    gsl_vector_set(funcDeriv, 1, secondDeriv);
}

void schrodingerEq_bound (double var, gsl_vector* funcVal, gsl_vector* funcDeriv){
    double thisFuncVal = gsl_vector_get(funcVal, 0);
    double firstDeriv = gsl_vector_get(funcVal, 1);
    double secondDeriv = (-2)*(1.0/var + bound_energy)*thisFuncVal;

    gsl_vector_set(funcDeriv, 0, firstDeriv);
    gsl_vector_set(funcDeriv, 1, secondDeriv);
}

void wavefunc(gsl_vector* vals, gsl_vector* funcVals){
    int dim  =  2;                                          // The order of the harmonic function
    gsl_vector* funcValRight  =  gsl_vector_alloc(dim);     // Vector to hold final value
    gsl_vector* funcValLeft   =  gsl_vector_calloc(dim);    // Vector to hold initial value

    // The ODE is defined on the interval [ leftEndpt ,  rightEndpt ]
    double  leftEndpt   =   1e-3;                               // ODE diverges at origin, so we choose small values
    double  rightEndpt  =   maxPt;
    double  absAcc      =   1e-3;                               // Absolute accuracy
    double  relAcc      =   1e-3;                               // Relative accuracy
    double  step        =   (rightEndpt - leftEndpt) / 10;      // Initial stepsize
    energy              =   gsl_vector_get(vals, 0);

    gsl_vector_set(funcValLeft, 0, (leftEndpt - leftEndpt*leftEndpt));
    gsl_vector_set(funcValLeft, 1, (1 - 2*leftEndpt));

    rkdriver(schrodingerEq, leftEndpt, funcValLeft, rightEndpt, funcValRight, step, absAcc, relAcc, NULL);
    gsl_vector_set(funcVals, 0, gsl_vector_get(funcValRight,0) );
}

void wavefunc_bound(gsl_vector* vals, gsl_vector* funcVals){
    int dim  =  2;                                          // The order of the harmonic function
    gsl_vector* funcValRight  =  gsl_vector_alloc(dim);     // Vector to hold final value
    gsl_vector* funcValLeft   =  gsl_vector_calloc(dim);    // Vector to hold initial value

    // The ODE is defined on the interval [ leftEndpt ,  rightEndpt ]
    double  leftEndpt   =   1e-5;
    double  rightEndpt  =   bound_maxPt;
    double  absAcc      =   1e-5;                               // Absolute accuracy
    double  relAcc      =   1e-5;                               // Relative accuracy
    double  step        =   (rightEndpt - leftEndpt) / 10;      // Initial stepsize
    bound_energy        =   gsl_vector_get(vals, 0);

    gsl_vector_set(funcValLeft, 0, (leftEndpt - leftEndpt*leftEndpt));
    gsl_vector_set(funcValLeft, 1, (1 - 2.*leftEndpt));

    rkdriver(schrodingerEq_bound, leftEndpt, funcValLeft, rightEndpt, funcValRight, step, absAcc, relAcc, NULL);
    gsl_vector_set(funcVals, 0, gsl_vector_get(funcValRight, 0) - rightEndpt*exp(-sqrt((-2)*bound_energy)*rightEndpt));
}

int main(int argc, char* argv[]){
    printf( "\nPart A) ");
    printf("\nRoot finding method on  f(x) = 1 + (x-a)² + (y-b)² ...\n");
    printf("Using a = %i, b = %i\n", 6, 13);

    int numOfDims = 2;
    double tolerance = 1e-5;
    gsl_vector* minimum = gsl_vector_alloc(numOfDims);
    gsl_vector_set(minimum, 0, 4);
    gsl_vector_set(minimum, 1, 10);

    printf("Tolerance: E = %g\n", tolerance);
    printf("Initial guess: (x, y)  = (%g, %g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));
    printf("Minimum:(x, y) = (6, 13)\n");
    newtonMethod(testFunc, minimum, tolerance);
    printf("Found minima: (x, y)  = (%g, %g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));


    printf("\nRoot finding on Rosenbrock Valley function\n");

    gsl_vector_set(minimum, 0, 0.5);
    gsl_vector_set(minimum, 1, 0.5);

    printf("Using tolerance: E  = %g\n", tolerance);
    printf("Initial guess: (x, y)  = (%g, %g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));
    printf("Minimum: (x, y)  = (1, 1)\n");
    newtonMethod(rosenbrockValley_grad, minimum, tolerance);
    printf("Found minima: (x, y)  = (%g, %g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));


    printf( "\nPart B)" );

    gsl_vector* minimum_hydrogen        =   gsl_vector_alloc(1);
    gsl_vector* minimum_hydrogen_bound  =   gsl_vector_alloc(1);
    gsl_vector_set(minimum_hydrogen,       0, -3);
    gsl_vector_set(minimum_hydrogen_bound, 0, -1);

    maxPt = 8.0;
    bound_maxPt = 0.5;
    newtonMethod(wavefunc,       minimum_hydrogen,       tolerance);
    newtonMethod(wavefunc_bound, minimum_hydrogen_bound, tolerance);

    printf("Energy = %g\n", energy);

    int dim  =  2;                                          // The order of the harmonic function
    gsl_vector* funcValRight  =  gsl_vector_alloc(dim);     // Vector to hold final value
    gsl_vector* funcValLeft   =  gsl_vector_calloc(dim);    // Vector to hold initial value

    // The ODE is defined on the interval [ leftEndpt ,  rightEndpt ]
    double  leftEndpt   =   1e-3;
    double  rightEndpt  =   maxPt;
    double  absAcc      =   1e-3;                               // Absolute accuracy
    double  relAcc      =   1e-3;                               // Relative accuracy
    double  step        =   (rightEndpt - leftEndpt) / 10;      // Initial stepsize

    gsl_vector_set(funcValLeft, 0, (leftEndpt - leftEndpt*leftEndpt));
    gsl_vector_set(funcValLeft, 1, (1 - 2.*leftEndpt));

    FILE* path2File = fopen(argv[1], "w"); // Set up filestream to write ODE solution to
    rkdriver(schrodingerEq, leftEndpt, funcValLeft, rightEndpt, funcValRight, step, absAcc, relAcc, path2File);
    fclose(path2File);

    printf( "\nPart C)");
    printf("Energy= %g\n", bound_energy);

    double exact_energy = -0.5;
    FILE* convergenceData = fopen("convergence.txt", "w");
    for (double pt = 0.1; pt <= 8.0; pt += 0.1) {
        maxPt = pt;
        bound_maxPt = pt;

        gsl_vector_set(minimum_hydrogen,       0, -3);
        gsl_vector_set(minimum_hydrogen_bound, 0, -1);
        newtonMethod(wavefunc,       minimum_hydrogen,       tolerance);
        newtonMethod(wavefunc_bound, minimum_hydrogen_bound, tolerance);

        fprintf(convergenceData, "%g\t%g\t%g\n", maxPt, fabs(energy - exact_energy), fabs(bound_energy - exact_energy));
    }
    fclose(convergenceData);

    gsl_vector_free(minimum);
    gsl_vector_free(minimum_hydrogen);
    gsl_vector_free(minimum_hydrogen_bound);


    return 0;
}
