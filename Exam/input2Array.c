#include <stdio.h>
#include "input2Array.h"
#include <stdlib.h>

void input2Array(double* XData, double* YData, char* inputFilename ){
  int numOfReturnVals = 2;  // The return value of fscanf() is the number of fetched variables

  int    input = numOfReturnVals; 			// Input will receive the return value from fscanf
  double argX;						              // Variable to hold the X data from stdin
  double argY;						              // Variable to hold the X data from stdin

  FILE* myInputFileStream  = fopen(inputFilename, "r");

  int id = 0;
  while( input != EOF) { // While we are not at the end of the input stream
    if ( input == numOfReturnVals){		 // If input has returned success
      input = fscanf(myInputFileStream, "%lg\t%lg", &argX, &argY); // Fetch data from stream and place it at the address of argX

      XData[id] = argX;
      YData[id] = argY;

      id++;
    }
    else { // If we are not successfull, say if there was an error with fscanf()
      fprintf(stderr, "Failed to read input.\n"); // Print this to stderr so we know
      exit(-1); 																	// Terminate program
    }
  }

  fclose(myInputFileStream);
}
