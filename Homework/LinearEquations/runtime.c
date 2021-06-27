#include <time.h>
#include <stdio.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include "timer.h"
#include "gramSchmidt.h"
#include "extraFuncs.h"


void runtime(int Reps, int startRep, char* outputFilename, char* gsl_outputFilename, unsigned int* seed){
  double scale      =  0   ;
  int Dims;

  FILE* myOutputFileStream      =  fopen(outputFilename,  "w");
  FILE* myOutputFileStream_gsl  =  fopen(gsl_outputFilename, "w");


  for (int rep = startRep; rep < Reps; rep++){
    Dims = rep;

    gsl_matrix* my_ortg    =  gsl_matrix_alloc(Dims, Dims);
    gsl_matrix* my_triang  =  gsl_matrix_alloc(Dims, Dims);
    gsl_vector* my_vecTmp  =  gsl_vector_alloc(Dims           );

    set_data(my_ortg, my_vecTmp, seed);

    clock_t my_begin  = clock(); 
    clock_t my_end    = clock(); 

    my_begin = clock(); 

    
    gramSchmidt_decomp(my_ortg, my_triang);

    my_end = clock(); 

    if (rep == startRep){
      scale = (double)(timer(my_end, my_begin));
    }

    fprintf(myOutputFileStream, "%d\t%g\t%g\n", Dims, (double)(timer(my_end, my_begin)), pow(((double)Dims)/startRep, 3)*scale);

    gsl_matrix_free(my_ortg);
    gsl_matrix_free(my_triang);
    gsl_vector_free(my_vecTmp);
  }

  for (int rep = startRep; rep < Reps; rep++){
    Dims = rep;

    gsl_matrix* gsl_ortg    =  gsl_matrix_alloc(Dims, Dims);
    gsl_matrix* gsl_triang  =  gsl_matrix_alloc(Dims, Dims);
    gsl_vector* my_vecTmp   =  gsl_vector_alloc(Dims);

    set_data(gsl_ortg, my_vecTmp , seed);

    clock_t gsl_begin = clock();
    clock_t gsl_end   = clock();

    gsl_begin = clock();
    gsl_linalg_QR_decomp(gsl_ortg, my_vecTmp);
    gsl_end = clock();

    fprintf(myOutputFileStream_gsl, "%d\t%g\n", Dims, (double)(timer(gsl_end, gsl_begin)));

    gsl_matrix_free(gsl_ortg);
    gsl_matrix_free(gsl_triang);
    gsl_vector_free(my_vecTmp);
  }
  fclose(myOutputFileStream);
  fclose(myOutputFileStream_gsl);
}
