//
// Created by marc on 5/6/21.
//


#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<float.h>

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.22045e-16
#endif

#include "minimization.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void numeric_gradient (double func(gsl_vector*), gsl_vector* minimum, gsl_vector* gradient){
    double stepSize =   sqrt(DBL_EPSILON);

    double funcVal  =   func(minimum);
    int numOfDims   =   minimum -> size;

    for(int dimId=0; dimId < numOfDims; ++dimId){
        double step;
        double minimum_i   =   gsl_vector_get(minimum, dimId);

        if (fabs(minimum_i) < stepSize) {
            step = stepSize;
        }
        else {
            step = fabs(minimum_i) * stepSize;
        }

        gsl_vector_set(minimum,  dimId,  minimum_i + step                 );
        gsl_vector_set(gradient, dimId,  (func(minimum) - funcVal) / step );
        gsl_vector_set(minimum,  dimId,  minimum_i - step                 );
    }
}

void quasi_newtonMethod(double func(gsl_vector*), gsl_vector* minimum, double tolerance){
    double stepSize = sqrt(DBL_EPSILON);
    int dimensions = minimum -> size;

    int numOfSteps      =   0;
    int numOfResets     =   0;
    int numOfScales     =   0;

    gsl_matrix* inverse_hessianMatrix   =   gsl_matrix_alloc(dimensions, dimensions);
    gsl_matrix_set_identity(inverse_hessianMatrix);

    gsl_matrix* identity   =   gsl_matrix_alloc(dimensions, dimensions);
    gsl_matrix_set_identity(identity);

    gsl_vector* gradientVal             =   gsl_vector_alloc(dimensions);
    gsl_vector* newtonStep              =   gsl_vector_alloc(dimensions);
    gsl_vector* minimum_next            =   gsl_vector_alloc(dimensions);
    gsl_vector* gradientVal_next        =   gsl_vector_alloc(dimensions);
    gsl_vector* solution                =   gsl_vector_alloc(dimensions);
    gsl_vector* solutionChange          =   gsl_vector_alloc(dimensions);
    gsl_vector* broydenVec              =   gsl_vector_alloc(dimensions);

    numeric_gradient(func, minimum, gradientVal);
    double funcVal  =  func(minimum);
    double funcVal_next;

    while(numOfSteps < 1000){
        numOfSteps++;

        gsl_blas_dgemv(CblasNoTrans, -1, inverse_hessianMatrix, gradientVal, 0, newtonStep);
        if( gsl_blas_dnrm2(newtonStep) < stepSize * gsl_blas_dnrm2(minimum) ) {
            fprintf(stderr,"quasi_newtonMethod: |Dx| < stepSize * |x|\n");
            break;
        }
        if( gsl_blas_dnrm2(gradientVal) < tolerance ) {
            fprintf(stderr,"quasi_newtonMethod: |grad| < acc\n");
            break;
        }

        double scale = 1;

        while(1){
            gsl_vector_memcpy(minimum_next, minimum);
            gsl_vector_add(minimum_next, newtonStep);

            funcVal_next = func(minimum_next);

            double sTransGrad;
            gsl_blas_ddot(newtonStep, gradientVal, &sTransGrad);

            if(funcVal_next < funcVal + 0.01 * sTransGrad){
                numOfScales++;
                break;
            }
            if(scale < stepSize){
                numOfResets++;
                gsl_matrix_set_identity(inverse_hessianMatrix);
                break;
            }
            scale*=0.5;
            gsl_vector_scale(newtonStep, 0.5);
        }

        numeric_gradient(func, minimum_next, gradientVal_next); // Compute gradient in the next step

        gsl_vector_memcpy(solution, gradientVal_next);
        gsl_blas_daxpy(-1, gradientVal, solution); /* y=grad(x+s)-grad(x) */
        gsl_vector_memcpy(solutionChange, newtonStep); /* u=s */
        gsl_blas_dgemv(CblasNoTrans, -1, inverse_hessianMatrix, solution, 1, solutionChange); /* u=s-By */

        gsl_matrix* solChangeSolChangeTrans = gsl_matrix_calloc(dimensions, dimensions); //u*u^T
        gsl_blas_dsyr(CblasUpper, 1.0, solutionChange, solChangeSolChangeTrans);
        double solChangeTransSol; // u^T*y
        gsl_blas_ddot(solutionChange, solution, &solChangeTransSol);
        if(fabs(solChangeTransSol) > 1e-12){ // SR1 update
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0 / solChangeTransSol, solChangeSolChangeTrans, identity, 1.0, inverse_hessianMatrix); // B= B + delta B
        }

        gsl_vector_memcpy(minimum, minimum_next);
        gsl_vector_memcpy(gradientVal, gradientVal_next);
        funcVal = funcVal_next;
    }

    gsl_matrix_free(inverse_hessianMatrix);
    gsl_matrix_free(identity);
    gsl_vector_free(gradientVal);
    gsl_vector_free(newtonStep);
    gsl_vector_free(minimum_next);
    gsl_vector_free(gradientVal_next);
    gsl_vector_free(solution);
    gsl_vector_free(solutionChange);
    gsl_vector_free(broydenVec);

    fprintf(stderr, "quasi_newtonMethod: \n  number of steps = %i,\n  number of scales = %i,\n  number of hessian matrix resets = %i\n  f(x) = %.1e\n\n", numOfSteps, numOfScales, numOfResets, funcVal);

}