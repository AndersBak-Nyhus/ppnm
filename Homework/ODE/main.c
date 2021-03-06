#include <math.h>
#include <gsl/gsl_matrix.h>
#include "RK-ODE.h"
#include "functions.h"

int main(int argc, char* argv[]){
    if (argc < 2){  // Check that we passed any arguments
        fprintf(stderr, "No arguments passed.");
        exit(-1);
    }


    int harmonicDim  =  2;                                                  
    
    gsl_vector* harmonicFuncValRight  =  gsl_vector_alloc(harmonicDim); 
    
    gsl_vector* harmonicFuncValLeft   =  gsl_vector_calloc(harmonicDim);   
    
    gsl_vector_set(harmonicFuncValLeft, 1, 1);                         
    
    
    double  leftEndpt   =   0.0;
    double  rightEndpt  =   2*M_PI;
    double  absAcc      =   1e-3;                               
    double  relAcc      =   1e-3;                               
    double  step        =   (rightEndpt - leftEndpt) / 10;      

    FILE* harmonic_outputStream = fopen(argv[1], "w"); 
    
    rkdriver(&harmonicFunc, leftEndpt, harmonicFuncValLeft, rightEndpt, harmonicFuncValRight, step, absAcc, relAcc, harmonic_outputStream);    
    fclose(harmonic_outputStream);
     
     

    leftEndpt   =   0.0;
    rightEndpt  =   100;
    int SIRdim  =   3;
    gsl_vector* SIRFuncValRight  =  gsl_vector_alloc(  SIRdim );
    gsl_vector* SIRFuncValLeft   =  gsl_vector_calloc( SIRdim );


	//part B)
    //covid-19 data
    double populationSize   =   5808180;
    double wasInfected      =   237792;
    double recovered        =   226630;
    double isInfected       =   wasInfected - recovered;
    double dead             =   2441;
    double vaccinated       =   445566;
    double removed          =   dead + recovered + vaccinated;

    
    gsl_vector_set(SIRFuncValLeft, 0, populationSize - isInfected - removed);
    gsl_vector_set(SIRFuncValLeft, 1, isInfected);
    gsl_vector_set(SIRFuncValLeft, 2, removed);

    FILE* SIR_outputStream = fopen(argv[2], "w");
    rkdriver(&SIRmodel, leftEndpt, SIRFuncValLeft, rightEndpt, SIRFuncValRight, step, absAcc, relAcc, SIR_outputStream);    
    fclose(SIR_outputStream);

    gsl_vector* SIR2funcValRight  =  gsl_vector_alloc(SIRdim);
    gsl_vector* SIR2funcValLeft   =  gsl_vector_calloc(SIRdim);
    gsl_vector_set(SIR2funcValLeft, 0, populationSize - isInfected - removed);
    gsl_vector_set(SIR2funcValLeft, 1, isInfected);
    gsl_vector_set(SIR2funcValLeft, 2, removed);
    FILE* SIR2_outputStream = fopen(argv[3], "w");
    rkdriver(&SIRmodel2, leftEndpt, SIR2funcValLeft, rightEndpt, SIR2funcValRight, step, absAcc, relAcc, SIR2_outputStream);    
    fclose(SIR2_outputStream);


	//part C)
    
    leftEndpt           =   0.0;
    rightEndpt          =   6.32591398;
    int threebodyDim    =   12;

    gsl_vector* threebodyFuncValRight  =  gsl_vector_alloc(threebodyDim);
    gsl_vector* threebodyFuncValLeft   =  gsl_vector_calloc(threebodyDim);

    double init_pos_x_1 =  0.97000436;
    double init_pos_y_1 = -0.24308753;
    double init_pos_x_2 = -0.97000436;
    double init_pos_y_2 =  0.24308753;
    double init_pos_x_3 =  0 ;
    double init_pos_y_3 =  0 ;
    double init_vel_x_3 = -0.93240737;
    double init_vel_y_3 = -0.86473146;
    double init_vel_x_1 = -init_vel_x_3/2;
    double init_vel_y_1 = -init_vel_y_3/2;
    double init_vel_x_2 = -init_vel_x_3/2;
    double init_vel_y_2 = -init_vel_y_3/2;

    gsl_vector_set(threebodyFuncValLeft, 0, init_pos_x_1);
    gsl_vector_set(threebodyFuncValLeft, 1, init_pos_y_1);
    gsl_vector_set(threebodyFuncValLeft, 2, init_pos_x_2);
    gsl_vector_set(threebodyFuncValLeft, 3, init_pos_y_2);
    gsl_vector_set(threebodyFuncValLeft, 4, init_pos_x_3);
    gsl_vector_set(threebodyFuncValLeft, 5, init_pos_y_3);
    gsl_vector_set(threebodyFuncValLeft, 6, init_vel_x_1);
    gsl_vector_set(threebodyFuncValLeft, 7, init_vel_y_1);
    gsl_vector_set(threebodyFuncValLeft, 8, init_vel_x_2);
    gsl_vector_set(threebodyFuncValLeft, 9, init_vel_y_2);
    gsl_vector_set(threebodyFuncValLeft, 10, init_vel_x_3);
    gsl_vector_set(threebodyFuncValLeft, 11, init_vel_y_3);

    FILE* threebody_outputStream = fopen(argv[4], "w");
    rkdriver(&threeBodyProb, leftEndpt, threebodyFuncValLeft, rightEndpt, threebodyFuncValRight, step, absAcc, relAcc, threebody_outputStream);   
    fclose(threebody_outputStream);


    // Free allocated memory
    gsl_vector_free(harmonicFuncValRight);
    gsl_vector_free(harmonicFuncValLeft);
    gsl_vector_free(SIRFuncValRight);
    gsl_vector_free(SIRFuncValLeft);
    gsl_vector_free(SIR2funcValRight);
    gsl_vector_free(SIR2funcValLeft);
    gsl_vector_free(threebodyFuncValRight);
    gsl_vector_free(threebodyFuncValLeft);

    return 0;
}

