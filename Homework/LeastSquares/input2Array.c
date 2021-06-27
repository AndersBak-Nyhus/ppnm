#include <stdio.h>
#include <stdlib.h>
#include "input2Array.h"

void input2Array(double* x, double* y, char* InputFile ){
  int RetVals = 3;  

  int    input = RetVals; 			
  
  double argX;						              
  
  double argY;
  
  double argErr;

  FILE* myInputFileStream  = fopen(InputFile, "r");

  int id = 0;
  while( input != EOF) { 
    if ( input == RetVals){
      input = fscanf(myInputFileStream, "%lg\t%lg\t%lg", &argX, &argY, &argErr); 

      x[id] = argX;
      y[id] = argY;

      id++;
    }
    else { 
      fprintf(stderr, "Failed to read input.\n");
      exit(-1); 																	// Terminate program
    }
  }

  fclose(myInputFileStream);
}
