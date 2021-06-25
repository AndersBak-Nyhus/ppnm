#include <float.h>
#include <assert.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "rootfinding.h"
#include "gramSchmidt.h"

void newtonMethod(void func(gsl_vector* point, gsl_vector* funcVals), gsl_vector* startingPoint, double tolerance){
    double stepSize = sqrt(DBL_EPSILON);
    int dimensions = startingPoint -> size;
    int count = 0;

    gsl_matrix* jacobianMatrix      =   gsl_matrix_alloc(dimensions, dimensions );
    gsl_vector* funcVal             =   gsl_vector_alloc(dimensions             );
    gsl_vector* funcVal_tmp         =   gsl_vector_alloc(dimensions             );
    gsl_matrix* triangMat           =   gsl_matrix_alloc(dimensions, dimensions );
    gsl_vector* solution            =   gsl_vector_alloc(dimensions             );
    gsl_vector* solution_scaled     =   gsl_vector_alloc(dimensions             );
    gsl_vector* nextPoint           =   gsl_vector_alloc(dimensions             );
    gsl_vector* nextFuncVal         =   gsl_vector_alloc(dimensions             );

    func(startingPoint, funcVal);
    while (gsl_blas_dnrm2(funcVal) > tolerance) {
        count++;
        assert(count < 1e4);

        for (int dimId = 0; dimId < dimensions; dimId++){
            gsl_vector_set(startingPoint, dimId, gsl_vector_get(startingPoint, dimId) + stepSize);
            func(startingPoint, funcVal_tmp);

            for (int subDimId = 0; subDimId < dimensions; subDimId++){

                double funcVals_tmp_id  =   gsl_vector_get(funcVal_tmp, subDimId);
                double funcVals_id      =   gsl_vector_get(funcVal, subDimId);
                double funcValDiff      =   funcVals_tmp_id - funcVals_id;
                double matrixElement    =   funcValDiff / stepSize;

                gsl_matrix_set(jacobianMatrix, subDimId, dimId, matrixElement);
            }
            gsl_vector_set(startingPoint, dimId, gsl_vector_get(startingPoint, dimId) - stepSize);
        }
        gsl_vector_scale(funcVal, -1.0);
        gramSchmidt_decomp(jacobianMatrix, triangMat);
        gramSchmidt_solve(jacobianMatrix, triangMat, funcVal, solution);
        gsl_vector_scale(funcVal, -1.0);

        double scale = 2;
        while ((gsl_blas_dnrm2(nextFuncVal) >= (1 - scale / 2) * gsl_blas_dnrm2(funcVal) ) && scale >= 0.02 ){
            scale /= 2;
            gsl_vector_memcpy(solution_scaled, solution);
            gsl_vector_scale(solution_scaled, scale);
            gsl_vector_memcpy(nextPoint, solution_scaled);
            gsl_vector_add(nextPoint, startingPoint);

            func(nextPoint, nextFuncVal);
        }
        gsl_vector_memcpy(startingPoint, nextPoint);
        gsl_vector_memcpy(funcVal, nextFuncVal);

        if (gsl_blas_dnrm2(solution) < stepSize){
            break;
        }
    }

    gsl_vector_free(funcVal);
    gsl_vector_free(nextFuncVal);
    gsl_vector_free(funcVal_tmp);
    gsl_vector_free(solution);
    gsl_vector_free(solution_scaled);
    gsl_vector_free(nextPoint);
    gsl_matrix_free(triangMat);
    gsl_matrix_free(jacobianMatrix);
}
