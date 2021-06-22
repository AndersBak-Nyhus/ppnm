#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "extraFuncs.h"
//creates random number
double randomNumber( unsigned int *seed ){
  double maxRand      =   (double)RAND_MAX;           
  double randNum      =   (double)rand_r( seed );     
  return randNum/maxRand;
}
//print vector in output
void vector_print(char* string, gsl_vector* vector){
	printf("%s\n", string);
	for(int iter = 0; iter < vector -> size; iter++){
		printf("%10g ", gsl_vector_get(vector, iter));
	}
  printf("\n");
}

//setting data inside matrices
void set_data(gsl_matrix* tallMatrix, gsl_vector* tallRHSVector, unsigned int *seed){

  for(int rowId = 0; rowId < tallMatrix -> size1; rowId++){
		for(int colId = 0; colId < tallMatrix -> size2; colId++){
  		gsl_matrix_set(tallMatrix, rowId, colId, randomNumber(seed));
		}
  	gsl_vector_set(tallRHSVector, rowId, randomNumber(seed));
  }
}

//printing matrix in output
void print_matrix(int Rows, gsl_matrix* matrixToPrint, char* string ){
  printf("\n%s\n", string);
  for (int rowId = 0; rowId < Rows; rowId++){
    gsl_vector_view matrixToPrint_row = gsl_matrix_row (matrixToPrint, rowId);
    gsl_vector* vector = &matrixToPrint_row.vector;
  	for(int iter = 0; iter < vector -> size; iter++){
      if ( gsl_vector_get(vector, iter) > 1e-10 ){
  		    printf("%10g\t", gsl_vector_get(vector, iter));
      }
      else { printf("%10g\t", 0.0); }
  	}
    printf("\n");
  }
}
