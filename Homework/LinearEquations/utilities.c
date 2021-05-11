#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "utilities.h"

double randomNumber( unsigned int *seed ){
  double maxRand      =   (double)RAND_MAX;           // Maximum random number, cast to double
  double randNum      =   (double)rand_r( seed );     // Generate pseudo-random number from seed, cast to double
  return randNum/maxRand;
}

void vector_print(char* string, gsl_vector* vector){
	printf("%s\n", string);
	for(int iter = 0; iter < vector -> size; iter++){
		printf("%10g ", gsl_vector_get(vector, iter));
	}
  printf("\n");
}

void set_data(gsl_matrix* tallMatrix, gsl_vector* tallRHSVector, unsigned int *seed){
  /* Method to set the data of a tall matrix using pseudorandom numbers.
   *
   *  ¤ gsl_matrix*    tallMatrix     : An arbitrary n x m matrix (n >= m)
   *  ¤ gsl_vector*    tallRHSVector  : A vector of dimension n x 1
   *  ¤ unsigned int*  seed           : A seed for the pseudorandom number generator
   *
   */

  for(int rowId = 0; rowId < tallMatrix -> size1; rowId++){
		for(int colId = 0; colId < tallMatrix -> size2; colId++){
  		gsl_matrix_set(tallMatrix, rowId, colId, randomNumber(seed));
		}
  	gsl_vector_set(tallRHSVector, rowId, randomNumber(seed));
  }
}
/*
void set_data_square(gsl_matrix* testMatSquare, gsl_vector* RHSvecSquare, unsigned int *LHSseed, unsigned int *RHSseed){
  for(int rowId = 0; rowId < testMatSquare -> size1; rowId++){
		for(int colId = 0; colId < testMatSquare -> size2; colId++){
  		gsl_matrix_set(testMatSquare, rowId, colId, randomNumber(LHSseed));
		}
  	gsl_vector_set(RHSvecSquare, rowId, randomNumber(RHSseed));
  }
}
*/
void print_matrix(int numOfRows, gsl_matrix* matrixToPrint, char* string ){
  printf("\n%s\n", string);
  for (int rowId = 0; rowId < numOfRows; rowId++){
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
