//
// Created by Marc on 4/7/21.
//

#include <math.h>

#include "RK-ODE.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


void rkstep12 ( void        (*func)(double, gsl_vector*, gsl_vector*),
                double      var                                      ,
                gsl_vector* funcVal                                  ,
                double      step                                     ,
                gsl_vector* funcStep                                 ,
                gsl_vector* err                                        ){
    /*
     *  A simple mid-point method runge-kutta stepper method, advancing
     *  the solution of an ordinary differential equation (ODE) by a step
     *  h. The method computes the step
     *
     *      y_{i + 1} = y_i + h*k
     *
     *  where the tangent, k, is given from
     *
     *      k_0     = f(x_0, y_0),
     *      k_{1/2} = f(x_0 + (1/2)*h, y_0 + (1/2)*h*k_0),
     *
     */
    int order = funcVal -> size; // This is the order of the differential equation

    gsl_vector* tangent_0   =  gsl_vector_alloc(order); // This is k_0
    gsl_vector* tangent_12  =  gsl_vector_alloc(order); // k_{1/2}
    gsl_vector* tmpFuncVal  =  gsl_vector_alloc(order); // and finally, a temporary variable

    // For some reason gsl BLAS methods (see rk45 method below) don't work so we use for-loops instead
    func(var, funcVal, tangent_0); // Call function to fill out tangent_0, creating the RHS of the ODE

    // Computing (y_0 + (1/2)*h*k_0))
    // id variables are temporary (scoped) variables that allow us to do the vector addition
    for (int id = 0; id < order; ++id ){
        double funcVal_id     =  gsl_vector_get(funcVal,   id);
        double tangent_0_id   =  gsl_vector_get(tangent_0, id);
        double tmpFuncVal_id  =  funcVal_id + 0.5*step*tangent_0_id;
        gsl_vector_set(tmpFuncVal, id, tmpFuncVal_id);
    }
    func(var + 0.5*step, tmpFuncVal, tangent_12); // Now we can fill k_{1/2} = f(x_0 + (1/2)*h, y_0 + (1/2)*h*k_0)

    // Now we can advance the solution (y_{i + 1} = y_i + h*k)
    for (int id = 0; id < order; ++id ){
        double funcVal_id      =  gsl_vector_get(funcVal,    id);
        double tangent_12_id   =  gsl_vector_get(tangent_12, id);
        double tmpFuncStep_id  =  funcVal_id + step*tangent_12_id;
        gsl_vector_set(funcStep, id, tmpFuncStep_id );
    }

    // Compute error estimate, following dmitris method (in the example) from the book
    for (int id = 0; id < order; ++id ){
        double tangent_0_id   =  gsl_vector_get(tangent_0, id);
        double tangent_12_id  =  gsl_vector_get(tangent_12, id);
        double tmpErr_id      =  (tangent_0_id - tangent_12_id) * step / 2;
        gsl_vector_set(err, id, tmpErr_id);
    }

    // Free allocated memory again
    gsl_vector_free(tangent_0);
    gsl_vector_free(tangent_12);
    gsl_vector_free(tmpFuncVal);
}

void rkstep45(
        void (*func)(double, gsl_vector*, gsl_vector*),  /* the f from dy/dt = f(t,y) */
        double var,                                      /* the current value of the variable */
        gsl_vector* funcVal,                             /* the current value y(t) of the sought function */
        double step,                                     /* the step to be taken */
        gsl_vector* funcStep,                            /* output: y(t+h) */
        gsl_vector* err                                  /* output: error estimate */
){
    // THIS FUNCTION DOESN'T WORK PROPERLY!
    int order = funcVal -> size;

    gsl_vector* funcDeriv       =   gsl_vector_alloc(order);
    gsl_vector* tangent_0       =   gsl_vector_alloc(order);
    gsl_vector* tangent_1       =   gsl_vector_alloc(order);
    gsl_vector* tangent_2       =   gsl_vector_alloc(order);
    gsl_vector* tangent_3       =   gsl_vector_alloc(order);
    gsl_vector* tangent         =   gsl_vector_alloc(order);
    gsl_vector* tangent_tmp_1   =   gsl_vector_alloc(order);
    gsl_vector* tangent_tmp_2   =   gsl_vector_alloc(order);

    gsl_matrix* identity        =   gsl_matrix_alloc(order, order);
    gsl_matrix_set_identity(identity);

    func(var, funcVal, funcDeriv);
    gsl_vector_memcpy(tangent_0, funcDeriv);
    gsl_vector_memcpy(tangent_1, tangent_0);
    gsl_blas_dgemv(CblasNoTrans, 1.0, identity, funcVal, 0.5*step, tangent_1);
    func(var + 0.5*step, tangent_1, tangent_1);
    gsl_vector_memcpy(tangent_2, tangent_1);
    gsl_blas_dgemv(CblasNoTrans, 1.0, identity, funcVal, 0.5*step, tangent_2);
    func(var + 0.5*step, tangent_2, tangent_2);
    gsl_vector_memcpy(tangent_3, tangent_2);
    gsl_blas_dgemv(CblasNoTrans, 1.0, identity, funcVal, step, tangent_3);
    func(var + step, tangent_3, tangent_3);

    gsl_vector_memcpy(tangent_tmp_1, tangent_1);
    gsl_vector_memcpy(tangent_tmp_2, tangent_3);
    gsl_blas_dgemv(CblasNoTrans, 1.0/6.0, identity, tangent_0, 1.0/3.0, tangent_tmp_1);
    gsl_blas_dgemv(CblasNoTrans, 1.0/3.0, identity, tangent_2, 1.0/6.0, tangent_tmp_2);
    gsl_vector_memcpy(tangent, tangent_tmp_2);
    gsl_blas_dgemv(CblasNoTrans, 1.0, identity, tangent_tmp_1, 1.0, tangent);


    gsl_blas_dgemv(CblasNoTrans, 1.0, identity, funcVal, step, tangent);
    gsl_vector_memcpy(funcStep, tangent);

    gsl_vector_memcpy(err, tangent);
    gsl_blas_dgemv(CblasNoTrans, 1.0, identity, funcVal, -1.0, err);

    gsl_vector_free(tangent_0);
    gsl_vector_free(tangent_1);
    gsl_vector_free(tangent_2);
    gsl_vector_free(tangent_3);
    gsl_vector_free(tangent_tmp_1);
    gsl_vector_free(tangent_tmp_2);
    gsl_vector_free(tangent);
    gsl_vector_free(funcDeriv);
    gsl_matrix_free(identity);
}


void rkdriver(  void (*func)(double, gsl_vector*, gsl_vector*),
                double leftEndpt ,  gsl_vector* funcValLeft   ,
                double rightEndpt,  gsl_vector* funcValRight  ,
                double step      ,
                double absAcc    ,  double relAcc             ,
                FILE* path2File                                 ) {
    /*
     *  A Runge-Kutta driver method used to compute the solution to an
     *  ordinary differential equation (ODE) of the planar form
     *
     *      dy/dt  =  f(t, y)
     *
     *  Where y is a vector of dimensionality equal to the order. The ODE is
     *  defined on the interval [a, b]. This driver calls a stepper function
     *  e.g. rkstep45() or rkstep12() defined above, which advances the
     *  solution by a stepsize, h.
     *
     *      ¤  (*func)      :   The right hand side of the ODE, a void function
     *                          of signature (double, gsl_vector*, gsl_vector*)
     *      ¤  leftEndPt    :   The initial value, the left end point, a, of
     *                          the solution interval.
     *      ¤  funcValLeft  :   The initial function value y(a), the function
     *                          value at the left end of the solution interval.
     *      ¤  rightEndPt   :   The left end point, b, of the solution interval.
     *      ¤  funcValRight :   The computed function value y(b), the function
     *                          value at the right end of the solution interval.
     *                          A return value that is to be computed
     *      ¤  step         :   The initial step size h.
     *      ¤  absAcc       :   Absolute accuracy.
     *      ¤  relAcc       :   Relative accuracy
     *      ¤  path2File    :   A FILE* to a file stream, to which the driver
     *                          writes the path of the computed solution.
     *
     */

    int order = funcValLeft -> size;  // This is the order of the differential equation

    double  err         ;   // Variable to hold the error estimated from the norm of the error vector
    double  normFuncVal ;   // Variable to hold the norm of the function values vector
    double  tol         ;

    gsl_vector* funcStep     =  gsl_vector_alloc(order);    // Function value, y-vector, at next step
    gsl_vector* funcErr      =  gsl_vector_alloc(order);    // Function value error, dy-vector
    gsl_vector* thisFuncVal  =  gsl_vector_alloc(order);    // Function value, y-vector, at current step
    gsl_vector_memcpy(thisFuncVal, funcValLeft);            // Initiallize current y-vector as the intial value

    double pos = leftEndpt;
    while(pos < rightEndpt){ // Run through entire interval
        if (path2File != NULL){ // Write data to file (t, y, dy/dt)
            fprintf(path2File, "%.5g\t", pos);
            for (int id = 0; id < order; ++id){
                fprintf(path2File,"%.5g\t", gsl_vector_get(thisFuncVal, id));
            }
            fprintf(path2File, "%g\n", pos*exp(-pos));
        }

        double trueStep;        // The final stepsize to advance by
        double nextStep = step; // Start from initial step size passed as argument

        if (pos + nextStep > rightEndpt) {  // If we overshoot interval bounds
            nextStep = rightEndpt - pos;    // then define stepsize to be what remains
        }
        do {
            rkstep12(func, pos, thisFuncVal, nextStep, funcStep, funcErr);  // Advance with step

            err          =  gsl_blas_dnrm2(funcErr);    // Compute error
            normFuncVal  =  gsl_blas_dnrm2(funcStep);

            tol = (normFuncVal * relAcc + absAcc) * sqrt(nextStep / (rightEndpt - leftEndpt)); // Compute tolerance
            trueStep = nextStep;
            nextStep *= pow(tol / err, 0.25) * 0.95;    // Nextstep will be made smaller, and used in the next iteration if we did not converge

        } while(err > tol); // Keep testing for convergence within tolerance, otherwise keep making stepsize smaller until convergence

        gsl_vector_memcpy(thisFuncVal, funcStep);   // Next step will be the next current function value
        pos += trueStep;                            // Advance current position
    }

    gsl_vector_memcpy(funcValRight, funcStep);      // Finally fill out the function vlaue at the right end point.

    // Free allocated memory
    gsl_vector_free(funcStep);
    gsl_vector_free(funcErr);
    gsl_vector_free(thisFuncVal);
}
