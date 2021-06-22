#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "montecarlo.h"

void print_testResults(char* string, double integralVal, double exactVal, double integrationError){
    printf("\n%s  : %g\n", string, integralVal);
    printf("Actual error                     : %g\n", fabs(integralVal - exactVal));
    printf("Computed error estimate          : %g\n", integrationError);
}

void print_whichTest(char* string, double exactVal){
    printf("\n%s %g\n", string, exactVal);
    printf("-----------------------------------------------\n");
}

double debugFunc ( double* x ){
    return sqrt( x[0] ) ;
};

double testFunc ( double* x ){
    return 1/(M_PI*M_PI*M_PI*(1-cos(x[0])*cos(x[1])*cos(x[2])));
};


int main(int argc, char* argv[]){

    // ##########################################################
    // # ------------ TEST MONTE CARLO INTEGRATOR ------------- #
    // ##########################################################

    int numOfPts = (int)1e2;

    printf("-----------------------------------------------");
    printf("\nA)\nPlain Monte Carlo integration\n");
    printf("-----------------------------------------------\n");
    printf("Testing plain Monte Carlo integration routine...\n");

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
    print_testResults("The numerical estimate of the integral is", result_debug, exactVal, error_debug);

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

    print_testResults("The numerical estimate of the integral using pseudo-random numbers is", result_test, exactVal, error_test);

    printf("\n\n-----------------------------------------------");
    printf("\nB)\nQuasi random numbers\n");
    printf("-----------------------------------------------\n");
    printf("Testing quasi-random Monte Carlo integration routine...\n");
    print_testResults("Done! The numerical estimate of the integral using quasi-random  numbers is", result_test_quasi, exactVal, error_test_quasi);

    printf("Testing the scaling of errors in quasi- vs pseudo-random number sequences...\n");
    char filename[] = "error_scaling.txt";
    int totalNumOfReps  =  (int) 1e3;
    int stepSize        =  (int) 1e2;
    /*FILE* myOutPutFileStream = fopen(filename, "w");
    for (int rep = 1; rep <= totalNumOfReps; rep++){
        int thisNumOfPts = rep * stepSize;

        plain_montecarlo(       dim_test, lowerBound_test, upperBound_test, testFunc, thisNumOfPts, &result_test,       &error_test       );
        plain_montecarlo_quasi( dim_test, lowerBound_test, upperBound_test, testFunc, thisNumOfPts, &result_test_quasi, &error_test_quasi );

        fprintf(myOutPutFileStream, "%i\t%g\t%g\n", thisNumOfPts, error_test, error_test_quasi);
    }
    fclose(myOutPutFileStream);
*/    printf("Done! The data has been written to %s\n", filename);


    printf("\n\n-----------------------------------------------");
    printf("\nC)\nRecursive stratified sampling\n");
    printf("-----------------------------------------------\n");
    printf("Testing recursive stratified sampling Monte Carlo integration routine...\n");

    double absAcc = 1e-3;
    double relAcc = 1e-3;
    double error_test_strat = 0;

    int numOfRecalls = 0;
    double meanRecalls = 0;
    double result_test_strat = plain_montecarlo_stratifiedSampling( dim_test, testFunc, lowerBound_test, upperBound_test, absAcc, relAcc, numOfRecalls, meanRecalls );

    //  print_testResults("Done! The numerical estimate of the integral using quasi-random  numbers is", result_test_strat, exactVal, 0);


    printf("-----------------------------------------------\n");
    return 0;
}

