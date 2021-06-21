#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
int binsearch(int n, double* XData, double z){

	assert(n>1 && z>=XData[0] && z<=XData[n-1]);
	int i=0,j=n-1;

	while(j-i>1){
		int m = (i+j)/2;

		if(z>=XData[m]) i=m;
		else j=m;
	}
	return i;
}
