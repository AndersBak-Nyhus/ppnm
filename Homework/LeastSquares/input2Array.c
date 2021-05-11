#include <stdio.h>
#include <stdlib.h>
#include "input2Array.h"

void input2Array(double* x, double* y, char* InputFile ){
  int NumberOfReturnVals = 3;  // The return value of fscanf() is the number of fetched variables

  int    input = NumberOfReturnVals; 			// Input will receive the return value from fscanf
  double argX;						              // Variable to hold the X data from stdin
  double argY;						              // Variable to hold the X data from stdin
  double argErr;

  FILE* myInputFileStream  = fopen(InputFile, "r");

  int id = 0;
  while( input != EOF) { // While we are not at the end of the input stream
    if ( input == NumberOfReturnVals){		 // If input has returned success
      input = fscanf(myInputFileStream, "%lg\t%lg\t%lg", &argX, &argY, &argErr); // Fetch data from stream and place it at the address of argX

      x[id] = argX;
      y[id] = argY;

      id++;
    }
    else { // If we are not successfull, say if there was an error with fscanf()
      fprintf(stderr, "Failed to read input.\n"); // Print this to stderr so we know
      exit(-1); 																	// Terminate program
    }
  }

  fclose(myInputFileStream);
}
