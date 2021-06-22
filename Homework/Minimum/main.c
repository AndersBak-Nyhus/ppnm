#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_vector.h>

#include "minimization.h"

double energy;
double bound_energy;
double maxPt;
double bound_maxPt;

double testFunction(gsl_vector* vals){
    // 1 + (x-a)^2 + (y-b)^2
    double x = gsl_vector_get(vals, 0);
    double y = gsl_vector_get(vals, 1);

    double a = 6;
    double b = 13;

    double Vals = 1 + (x-a)*(x-a) + (y-b)*(y-b);

    return Vals;
}

double rosenbrockValley(gsl_vector* vals) {
    double x = gsl_vector_get(vals, 0);
    double y = gsl_vector_get(vals, 1);

    double Vals = (1-x)*(1-x)+100*(y-x*x)*(y-x*x);
    return Vals; 
}

double himmelblau(gsl_vector* vals){
	double x = gsl_vector_get(vals,0);
	double y = gsl_vector_get(vals,1);

	double Vals = (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
	return Vals;
}

int main(int argc, char* argv[]){

    printf("\nPart A)\n");

    printf("\nTest on f(x) = 1 + x² + y² ...\n");
    double min_x = 6;
    double min_y = 13;
    double ini_x = min_x - 2;
    double ini_y = min_y - 3;

    int dims = 2;
    double tolerance = 1e-5;
    gsl_vector* minimum = gsl_vector_alloc(dims);
    gsl_vector_set(minimum, 0, ini_x);
    gsl_vector_set(minimum, 0, ini_y);

    printf("Initial value (x, y) = (%g, %g)\n", ini_x, ini_y);
    quasi_newtonMethod(testFunction, minimum, tolerance);
    printf("Minima (x, y) = (%g, %g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));
    printf("Real minima (x, y) = (%g, %g)\n", min_x, min_y);


	printf("\n\nRosenbrock's\n");
	min_x = 1;
	min_y = 1;
	ini_x = min_x-1;
	ini_y = min_y-1;

	gsl_vector_set(minimum, 0, ini_x);
	gsl_vector_set(minimum, 0, ini_y);

    printf("Initial value (x, y) = (%g, %g)\n", ini_x, ini_y);
    quasi_newtonMethod(rosenbrockValley, minimum, tolerance);
    printf("Minima (x, y) = (%g, %g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));

        printf("\n\nHimmelblau's\n");
        min_x = 3;
        min_y = 2;
        ini_x = min_x-1;
        ini_y = min_y-1;

        gsl_vector_set(minimum, 0, ini_x);
        gsl_vector_set(minimum, 0, ini_y);

    printf("Initial value (x, y) = (%g, %g)\n", ini_x, ini_y);
    quasi_newtonMethod(himmelblau, minimum, tolerance);
    printf("Minima (x, y) = (%g, %g)\n", gsl_vector_get(minimum, 0), gsl_vector_get(minimum, 1));    
    return 0;
}
