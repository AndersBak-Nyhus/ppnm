#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "montecarlo.h"

void print_testResults(char* string, double integralVal, double exactVal, double integrationError){
    printf("\n%s  : %g\n", string, integralVal);
    printf("Actual error                     : %g\n", fabs(integralVal - exactVal));
    printf("Calculated error estimate          : %g\n", integrationError);
}

void print_whichTest(char* string, double exactVal){
    printf("\n%s %g\n", string, exactVal);
}

double debugFunc ( double* x ){
    return sqrt( x[0] ) ;
};

double testFunc ( double* x ){
    return 1/(M_PI*M_PI*M_PI*(1-cos(x[0])*cos(x[1])*cos(x[2])));
};


int main(int argc, char* argv[]){

    int numOfPts = (int)1e6;

    printf("\nPart A)\nMonte Carlo integration\n");
    printf("CheckMonte Carlo integration\n");

    double exactVal     =   2.0/3.0;
    print_whichTest("∫_0^1 dx √(x) = 2/3 =", exactVal);

    const int dim_debug = 1;
    double* lowerBound_debug = calloc( dim_debug, sizeof(double) );
    double* upperBound_debug = malloc( sizeof(double)*dim_debug );
    lowerBound_debug[0] = 0;
    upperBound_debug[0] = 1;
    double result_debug;
    double error_debug;

    result_debug = 0;
    error_debug = 0;

    plain_montecarlo(dim_debug, lowerBound_debug, upperBound_debug, debugFunc, numOfPts, &result_debug, &error_debug);
    print_testResults("Numerical estimate", result_debug, exactVal, error_debug);

    exactVal = 1.3932039296856768591842462603255;
    print_whichTest("∫_0^π dx/π ∫_0^π dy/π ∫_0^π  dz/π [1-cos(x)cos(y)cos(z)]^{-1} = Γ(1/4)4/(4π3) ≈", exactVal);

    const int dim_test = 3;
    double* lowerBound_test = calloc( dim_test, sizeof(double) );
    double* upperBound_test = malloc( sizeof(double)*dim_test );
    for (int axis = 0; axis < dim_test; axis++){
        lowerBound_test[axis] = 0;
        upperBound_test[axis] = M_PI;
    }
    double result_test;
    double error_test;
    double result_test_quasi;
    double error_test_quasi;

    result_test = 0;
    error_test = 0;
    result_test_quasi = 0;
    error_test_quasi = 0;

    plain_montecarlo(       dim_test, lowerBound_test, upperBound_test, testFunc, numOfPts, &result_test,       &error_test       );
    plain_montecarlo_quasi( dim_test, lowerBound_test, upperBound_test, testFunc, numOfPts, &result_test_quasi, &error_test_quasi );

    print_testResults("Numerical estimate with pseudo-random numbers", result_test, exactVal, error_test);

    printf("\nPart B)\nQuasi random numbers\n");
    printf("quasi-random Monte Carlo integration\n");
    print_testResults("Numerical estimate using quasi-random  numbers", result_test_quasi, exactVal, error_test_quasi);

    printf("errors in quasi- vs pseudo-random\n");
    char filename[] = "error.txt";
    int totalNumOfReps  =  (int) 1e3;
    int stepSize        =  (int) 1e2;




    printf("\nPart C)\nRecursive sampling\n");
    printf("Recursive sampling Monte Carlo integration\n");

    double absAcc = 1e-3;
    double relAcc = 1e-3;
    double error_test_strat = 0;

    int numOfRecalls = 0;
    double meanRecalls = 0;
    double result_test_strat = plain_montecarlo_stratifiedSampling( dim_test, testFunc, lowerBound_test, upperBound_test, absAcc, relAcc, numOfRecalls, meanRecalls );

    return 0;
}

