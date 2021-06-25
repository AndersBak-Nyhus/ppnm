#include <time.h>
#include <stdio.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include "timer.h"
#include "jacobi.h"
#include "extraFuncs.h"


void runtime(int Reps, int startRep, char* my_outputFilename, char* gsl_outputFilename, unsigned int* seed){
  double scale      =  0   ;
  int Dims;

  FILE* myOutputFileStream      =  fopen(my_outputFilename,  "w");
  FILE* myOutputFileStream_gsl  =  fopen(gsl_outputFilename, "w");


  for (int rep = startRep; rep < Reps; rep++){
    Dims = rep;

    gsl_matrix* symm    =  gsl_matrix_alloc(Dims, Dims);
    gsl_matrix* EigVec  =  gsl_matrix_alloc(Dims, Dims);
    
    
    set_data_symmetric(symm, seed);

    clock_t my_begin  = clock(); 
    clock_t my_end    = clock(); 

    my_begin = clock(); 

    
    jacobiDiag(symm, EigVec);

    my_end = clock(); 

    if (rep == startRep){
      scale = (double)(timer(my_end, my_begin));
    }

    fprintf(myOutputFileStream, "%d\t%g\t%g\n", Dims, (double)(timer(my_end, my_begin)), pow(((double)Dims)/startRep, 3)*scale);

    gsl_matrix_free(symm);
    gsl_matrix_free(EigVec);
  }

  for (int rep = startRep; rep < Reps; rep++){
    Dims = rep;

    gsl_matrix* gsl_symm    =  gsl_matrix_alloc(Dims, Dims);
    gsl_matrix* gsl_EigVec  =  gsl_matrix_alloc(Dims, Dims);
    gsl_vector* gsl_Diag   =  gsl_vector_alloc(Dims);

    set_data_symmetric(gsl_symm, seed);

    clock_t gsl_begin = clock();
    clock_t gsl_end   = clock();

    gsl_begin = clock();
    gsl_linalg_SV_decomp_jacobi(gsl_symm, gsl_EigVec, gsl_Diag);
    gsl_end = clock();

    fprintf(myOutputFileStream_gsl, "%d\t%g\n", Dims, (double)(timer(gsl_end, gsl_begin)));

    gsl_matrix_free(gsl_symm);
    gsl_matrix_free(gsl_EigVec);
    gsl_vector_free(gsl_Diag);
  }
  fclose(myOutputFileStream);
  fclose(myOutputFileStream_gsl);
}
