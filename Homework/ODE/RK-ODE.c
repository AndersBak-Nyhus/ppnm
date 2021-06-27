#include <math.h>
#include "RK-ODE.h"
#include "functions.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


void rkstep12 ( void        (*func)(double, gsl_vector*, gsl_vector*),
                double      var                                      ,
                gsl_vector* funcVal                                  ,
                double      step                                     ,
                gsl_vector* funcStep                                 ,
                gsl_vector* err                                        ){

    int order = funcVal -> size; 

    gsl_vector* tangent_0   =  gsl_vector_alloc(order);
    gsl_vector* tangent_12  =  gsl_vector_alloc(order);
    gsl_vector* tmpFuncVal  =  gsl_vector_alloc(order); 


    func(var, funcVal, tangent_0); 

  
    for (int id = 0; id < order; ++id ){
        double funcVal_id     =  gsl_vector_get(funcVal,   id);
        double tangent_0_id   =  gsl_vector_get(tangent_0, id);
        double tmpFuncVal_id  =  funcVal_id + 0.5*step*tangent_0_id;
        gsl_vector_set(tmpFuncVal, id, tmpFuncVal_id);
    }
    func(var + 0.5*step, tmpFuncVal, tangent_12); 
    
    for (int id = 0; id < order; ++id ){
        double funcVal_id      =  gsl_vector_get(funcVal,    id);
        double tangent_12_id   =  gsl_vector_get(tangent_12, id);
        double tmpFuncStep_id  =  funcVal_id + step*tangent_12_id;
        gsl_vector_set(funcStep, id, tmpFuncStep_id );
    }

   
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
        void (*func)(double, gsl_vector*, gsl_vector*),  
        double var,                                      
        gsl_vector* funcVal,                             
        double step,                                     
        gsl_vector* funcStep,                            
        gsl_vector* err                                  
){
    
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


    int order = funcValLeft -> size;  

    double  err         ;   
    double  normFuncVal ;   
    double  tol         ;

    gsl_vector* funcStep     =  gsl_vector_alloc(order);    
    gsl_vector* funcErr      =  gsl_vector_alloc(order);    
    gsl_vector* thisFuncVal  =  gsl_vector_alloc(order);    
    gsl_vector_memcpy(thisFuncVal, funcValLeft);            

    double pos = leftEndpt;
    while(pos < rightEndpt){ 
        if (path2File != NULL){ 
            fprintf(path2File, "%.5g\t", pos);
            for (int id = 0; id < order; ++id){
                fprintf(path2File,"%.5g\t", gsl_vector_get(thisFuncVal, id));
            }
            if ( func == harmonicFunc ){ 
                fprintf(path2File,"%.5g\n", sin(pos));
            }
            else {
                fprintf(path2File, "\n");
            }
        }

        double trueStep;        
        double nextStep = step; 

        if (pos + nextStep > rightEndpt) {  
            nextStep = rightEndpt - pos;    
        }
        do {
            rkstep12(func, pos, thisFuncVal, nextStep, funcStep, funcErr);  

            err          =  gsl_blas_dnrm2(funcErr);    
            normFuncVal  =  gsl_blas_dnrm2(funcStep);

            tol = (normFuncVal * relAcc + absAcc) * sqrt(nextStep / (rightEndpt - leftEndpt)); 
            trueStep = nextStep;
            nextStep *= pow(tol / err, 0.25) * 0.95;    

        } while(err > tol); 

        gsl_vector_memcpy(thisFuncVal, funcStep);   
        pos += trueStep;                            
    }

    gsl_vector_memcpy(funcValRight, funcStep);      

    // Free allocated memory
    gsl_vector_free(funcStep);
    gsl_vector_free(funcErr);
    gsl_vector_free(thisFuncVal);
}
