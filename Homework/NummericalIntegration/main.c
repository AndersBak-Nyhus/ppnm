#include <math.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>


#include "integration.h"

void print_testResults(char* string, double integralVal, double exactVal, double absAcc, double relAcc, double integrationError, int numOfCalls){
    printf("\n%s  : %g\n", string, integralVal);
    printf("Error goal                       : %.25g\n", absAcc + fabs( exactVal ) * relAcc);
    printf("Actual error                     : %.25g\n", fabs(integralVal - exactVal));
    printf("Computed error estimate          : %.25g\n", integrationError);
    printf("Function was called %i times...\n", numOfCalls);
}

void print_whichTest(char* string, double exactVal){
    printf("\n%s %g", string, exactVal);
    printf("\n---------------------------------");
}


int main(int argc, char* argv[]){

    // ###########################################################
    // # -------------- TEST ADAPTIVE INTEGRATOR --------------- #
    // ###########################################################

    double leftEndPt    =   0;
    double rightEndPt   =   1;
    double absAcc       =   1e-3;
    double relAcc       =   1e-3;

    int     numOfCalls          =   0;
    double  integrationError    =   0;

    printf("\nPart C) about error computation is included in parts A), B) as well...\n");
    printf("\n__________________________________________________________________");

    printf("\nA)\nRecursive adaptive integrator");
    printf("\nTesting integrator on different integrals...\n");

    double exactVal     =   2.0/3.0;
    print_whichTest("∫_0^1 dx √(x) = 2/3 =", exactVal);
    double firstTestFunc ( double x ){
        numOfCalls++;
        return sqrt( x ) ;
    };

    double integralVal  =   integrate( firstTestFunc, leftEndPt, rightEndPt, absAcc, relAcc, &integrationError );

    print_testResults("Result of numerical integration", integralVal, exactVal, absAcc, relAcc, integrationError, numOfCalls);



    double secondTestFunc ( double x ){
        numOfCalls++;
        return 4 * sqrt( 1 - x*x );
    };
    numOfCalls = 0;
    integrationError = 0;
    integralVal = integrate( secondTestFunc, leftEndPt, rightEndPt, absAcc, relAcc, &integrationError);
    exactVal = M_PI;
    print_whichTest("∫_0^1 dx 4√(1-x²) = π =", exactVal);
    print_testResults("Result of numerical integration", integralVal, exactVal, absAcc, relAcc, integrationError, numOfCalls);

    printf("\n__________________________________________________________________");
    // _________________________________________________________________________________________________________________


    // ######################################################################
    // # -- OPEN QUADRATURE WITH CLENSHAW-CURTIS VARIABLE TRANSFORMATION -- #
    // ######################################################################

    printf("\nB)\nOpen quadrature with Clenshaw–Curtis variable transformation\n");

    exactVal = 2.0;
    print_whichTest("∫_0^1 dx 1/√(x) = 2 = ", exactVal);

    double thirdTestFunc ( double x ){
        numOfCalls++;
        return 1/sqrt(x);
    };

    numOfCalls = 0;
    integrationError = 0;
    integralVal = open_quad(thirdTestFunc, leftEndPt, rightEndPt, absAcc, relAcc, &integrationError);
    print_testResults("Result of numerical integration with Clenshaw-Curtis", integralVal, exactVal, absAcc, relAcc, integrationError, numOfCalls);

    numOfCalls = 0;
    integrationError = 0;
    integralVal = adapt(thirdTestFunc, leftEndPt, rightEndPt, absAcc, relAcc, &integrationError);
    print_testResults("Result of numerical integration without Clenshaw-Curtis", integralVal, exactVal, absAcc, relAcc, integrationError, numOfCalls);

    exactVal = -4.0;
    print_whichTest("∫_0^1 dx ln(x)/√(x) = -4 = ", exactVal);

    double fourthTestFunc ( double x ){
        numOfCalls++;
        return log(x)/sqrt(x);
    };

    numOfCalls = 0;
    integrationError = 0;
    integralVal = open_quad(fourthTestFunc, leftEndPt, rightEndPt, absAcc, relAcc, &integrationError);
    print_testResults("Result of numerical integration with Clenshaw-Curtis", integralVal, exactVal, absAcc, relAcc, integrationError, numOfCalls);

    numOfCalls = 0;
    integrationError = 0;
    integralVal = adapt(fourthTestFunc, leftEndPt, rightEndPt, absAcc, relAcc, &integrationError);
    print_testResults("Result of numerical integration without Clenshaw-Curtis", integralVal, exactVal, absAcc, relAcc, integrationError, numOfCalls);

    exactVal = M_PI;
    print_whichTest("∫_0^1 dx 4√(1-x²) = π = ", exactVal);

    numOfCalls = 0;
    integrationError = 0;
    integralVal = open_quad(secondTestFunc, leftEndPt, rightEndPt, absAcc, relAcc, &integrationError);
    print_testResults("Result of numerical integration with Clenshaw-Curtis", integralVal, exactVal, absAcc, relAcc, integrationError, numOfCalls);

    numOfCalls = 0;
    integrationError = 0;
    integralVal = adapt(secondTestFunc, leftEndPt, rightEndPt, absAcc, relAcc, &integrationError);
    print_testResults("Result of numerical integration without Clenshaw-Curtis", integralVal, exactVal, absAcc, relAcc, integrationError, numOfCalls);

    double gslTestFunc( double x, void* params ){
        params = NULL;
        return secondTestFunc(x);
    }
    gsl_function my_gsl_test_func;
    my_gsl_test_func.function = &gslTestFunc;
    my_gsl_test_func.params = NULL;
    size_t limit   =   999 ;
    gsl_integration_cquad_workspace* workspace = gsl_integration_cquad_workspace_alloc( limit );
    double result;
    double absError;
    size_t numOfEvals;

    integrationError = 0;
    gsl_integration_cquad(&my_gsl_test_func, leftEndPt, rightEndPt, absAcc, relAcc, workspace, &result, &absError, &numOfEvals );

    print_testResults("Result of numerical integration with GSL Clenshaw-Curtis", result, exactVal, absAcc, relAcc, absError, (int)numOfEvals);
    workspace = NULL;
    gsl_integration_cquad_workspace_free(workspace);

    printf("__________________________________________________________________\n");
    
    double err = 0;
    // _________________________________________________________________________________________________________________


    // ######################################################################
    // # ------------------------ INFINITE LIMITS ------------------------- #
    // ######################################################################

    printf("\nC)\nInfinite limits\n");

    double fifthTestFunc( double x ){
        numOfCalls++;
        return exp(-x*x);
    }
    exactVal = sqrt(M_PI);
    print_whichTest("∫_-inf^inf dx exp(-x²) = √π =", exactVal);
    numOfCalls = 0;
    integrationError = 0;
    integralVal = integrate(fifthTestFunc, -INFINITY, INFINITY, absAcc, relAcc, &integrationError);
    print_testResults("Result of numerical integration", integralVal, exactVal, absAcc, relAcc, integrationError, numOfCalls);

    double gslTestFunc2( double x, void* params ){
        params = NULL;
        return fifthTestFunc(x);
    }
    gsl_function my_gsl_test_func2;
    my_gsl_test_func2.function = &gslTestFunc2;
    my_gsl_test_func2.params = NULL;
    gsl_integration_workspace* workspace2 = gsl_integration_workspace_alloc( limit );
    result = 0;
    absError = 0;
    numOfEvals = 0;
    gsl_integration_qagi(&my_gsl_test_func2, absAcc, relAcc, limit, workspace2, &result, &absError);
    print_testResults("Result of numerical GSL integration", result, exactVal, absAcc, relAcc, absError, (int)numOfEvals);
    printf("(QAGI-methods cannot return number of evaluations...)\n");



    double sixthTestFunc( double x ){
        numOfCalls++;
        return 1/(1 + x*x);
    }

    exactVal = M_PI/2;
    print_whichTest("∫_0^inf dx 1/(1+x²) = π/2 =", exactVal);
    numOfCalls = 0;
    integrationError = 0;
    integralVal = integrate( sixthTestFunc, 0, INFINITY, absAcc, relAcc, &integrationError);
    print_testResults("Result of numerical integration", integralVal, exactVal, absAcc, relAcc, integrationError, numOfCalls);

    double gslTestFunc3( double x, void* params ){
        params = NULL;
        return sixthTestFunc(x);
    }
    gsl_function my_gsl_test_func3;
    my_gsl_test_func3.function = &gslTestFunc3;
    my_gsl_test_func3.params = NULL;
    gsl_integration_workspace* workspace3 = gsl_integration_workspace_alloc( limit );
    result = 0;
    absError = 0;
    numOfEvals = 0;
    gsl_integration_qagiu(&my_gsl_test_func3, 0, absAcc, relAcc, limit, workspace3, &result, &absError);
    print_testResults("Result of numerical GSL integration", result, exactVal, absAcc, relAcc, absError, (int)numOfEvals);
    printf("(QAGI-methods cannot return number of evaluations...)\n");



    exactVal = M_PI/2;
    print_whichTest("∫-inf,0 dx  1/(1+x²) = π/2 =", exactVal);
    numOfCalls = 0;
    integrationError = 0;
    integralVal = integrate( sixthTestFunc, -INFINITY, 0, absAcc, relAcc, &integrationError);
    print_testResults("Result of numerical integration", integralVal, exactVal, absAcc, relAcc, integrationError, numOfCalls);


    gsl_function my_gsl_test_func4;
    my_gsl_test_func4.function = &gslTestFunc3;
    my_gsl_test_func4.params = NULL;
    gsl_integration_workspace* workspace4 = gsl_integration_workspace_alloc( limit );
    result = 0;
    absError = 0;
    numOfEvals = 0;
    gsl_integration_qagil(&my_gsl_test_func4, 0, absAcc, relAcc, limit, workspace4, &result, &absError);
    print_testResults("Result of numerical GSL integration", result, exactVal, absAcc, relAcc, absError, (int)numOfEvals);
    printf("(QAGI-methods cannot return number of evaluations...)\n");



    return 0;
}

