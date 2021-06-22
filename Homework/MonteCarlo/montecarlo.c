#include "montecarlo.h"

#include <math.h>
#include <stdlib.h>
#include <assert.h>

#define fracl(x) ((x) -  floorl(x))
#define PI 3.1415926535897932384626433832795028841971693993751L
#define RND ((double)rand()/RAND_MAX)

double van_der_corput_sequence(int id , int base){
    double corput_num = 0;
    double coprime_base = (double)1/base;
    while (id > 0){
        corput_num += (id % base ) * coprime_base ;
        id /= base ;
        coprime_base /= base;
    }
    return corput_num ;
}

void halton_sequence(int id, int dim, double *vec){
    int base[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53} ;
    int max_dim = sizeof (base) / sizeof(int);
    assert(dim <= max_dim) ;
    for (int axis = 0; axis < dim ; axis++){
        vec[axis] = van_der_corput_sequence(id + 1, base[axis]);
    }
}

void halton_sequence_secondary(int id, int dim, double *vec){
    int base[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67} ;
    int max_dim = sizeof (base) / sizeof(int);
    assert(dim <= max_dim) ;
    for (int axis = 0; axis < dim ; axis++){
        vec[axis] = van_der_corput_sequence(id + 1, base[axis]);
    }
}


void lattice(int d, double* x) {
    static int dim = 0;
    static int n = -1;
    static long double *alpha;

    if (d < 0) {
        dim = -d;
        n = 0;
        alpha = (long double *) realloc(alpha, dim * sizeof(long double));

        for (int i = 0; i < dim; i++) {
            alpha[i] = fracl(sqrtl(PI + i));
        }
    } else if (d > 0) {
        n++;
        assert(d == dim && n > 0);
        for (int i = 0; i < dim; i++) {
            x[i] = fracl(n * alpha[i]);
        }
    } else if (alpha != NULL) {
        free(alpha);
    }

}

void rand_num (int dim , const double *lowerBound , const double *upperBound , double *nums ) {
    for (int axis = 0; axis < dim; axis++) {
        nums[axis] = lowerBound[axis] + RND * (upperBound[axis] - lowerBound[axis]);
    }
}

void rand_num_halton_corput (int id, int dim , const double *lowerBound , const double *upperBound , double *vec ) {
    halton_sequence(id, dim, vec);
    for (int axis = 0; axis < dim; axis++) {
        vec[axis] = lowerBound[axis] + vec[axis] * (upperBound[axis] - lowerBound[axis]);
    }
}

void rand_num_halton_corput_secondary (int id, int dim , const double *lowerBound , const double *upperBound , double *vec ) {
    halton_sequence_secondary(id, dim, vec);
    for (int axis = 0; axis < dim; axis++) {
        vec[axis] = lowerBound[axis] + vec[axis] * (upperBound[axis] - lowerBound[axis]);
    }
}

void plain_montecarlo (int dim, double *lowerBound , double *upperBound, double func(double* nums) , int numOfPts, double* result , double* error ){

    double volume = 1;
    for (int axis = 0; axis < dim ; axis++) {
        volume*= upperBound[axis] - lowerBound[axis];
    }

    double sum = 0;
    double sum_squared = 0;
    double funcVal;
    double vec[dim];

    for (int id = 0; id < numOfPts; id ++){
        rand_num(dim, lowerBound, upperBound, vec) ;
        funcVal = func(vec) ;
        sum +=funcVal ;
        sum_squared += funcVal * funcVal;
    }

    double average = sum / numOfPts;
    double variance = sum_squared / numOfPts - average * average ;
    *result = average * volume;
    *error = sqrt(variance / numOfPts) * volume;
}

void plain_montecarlo_quasi (int dim, double *lowerBound , double *upperBound, double func(double* nums) , int numOfPts, double* result , double* error ){

    double volume = 1;
    for (int axis = 0; axis < dim ; axis++) {
        volume*= upperBound[axis] - lowerBound[axis];
    }

    double sum_first = 0;
    double sum_second = 0;
    double funcVal_first;
    double funcVal_second;
    double vec_first[dim];
    double vec_second[dim];

    for (int id = 0; id < floor(numOfPts/2); id++){
        rand_num_halton_corput(           id, dim, lowerBound, upperBound, vec_first  ) ;
        rand_num_halton_corput_secondary( id, dim, lowerBound, upperBound, vec_second ) ;
        funcVal_first   =  func(vec_first) ;
        funcVal_second  =  func(vec_second) ;

        if(!isinf(funcVal_first) && !isinf(funcVal_second)) {
            sum_first  += funcVal_first;
            sum_second += funcVal_second;
        }
    }

    double average = (sum_second + sum_first) / numOfPts;

    *result = average * volume;
    *error  = volume * fabs(sum_first - sum_second) / (numOfPts);
}


double plain_montecarlo_stratifiedSampling( int     dim                         ,
                                            double  func(double* vec)  ,
                                            double* lowerBound                  ,
                                            double* upperBound                  ,
                                            double  absAcc                      ,
                                            double  relAcc                      ,
                                            int     numOfRecalls                ,
                                            double  meanRecall                      ) {

    int numOfPts = 16 * dim ;

    double volume = 1;
    for(int axis = 0; axis < dim; axis++) {
        volume *= upperBound [ axis] - lowerBound [ axis ] ;
    }

    int numToTheLeft[dim];
    int numToTheRight[dim];

    double vec[dim];
    double averageLeft[dim];
    double averageRight[dim];
    double average = 0;

    for(int axis = 0; axis < dim; axis++) {
        averageLeft[axis]         = 0 ;
        averageRight[axis]        = 0 ;
        numToTheLeft[axis]     = 0 ;
        numToTheRight[axis]    = 0 ;
    }

    for (int id = 0; id < numOfPts; id++){
        rand_num(dim, lowerBound, upperBound, vec);
        double funcVal = func(vec);
        average += funcVal ;

        for(int axis = 0; axis < dim; axis++) {
            double midPt = (lowerBound[axis] + upperBound[axis]) / 2;
            if(vec[axis] > midPt){
                numToTheRight[axis]++;
                averageRight[axis] += funcVal;
            }
            else {
                numToTheLeft[axis]++;
                averageLeft[axis] += funcVal;
            }
        }
    }
    average /= numOfPts;
    for(int axis = 0; axis < dim; axis++) {
        averageLeft[axis]     /=  numToTheLeft[axis];
        averageRight[axis]    /=  numToTheRight[axis];
    }

    int axisDiv = 0;
    double maxVariance = 0;

    for(int axis = 0; axis < dim; axis++) {
        double variance = fabs(averageRight[axis] - averageLeft[axis]);
        if(variance > maxVariance) {
            maxVariance = variance;
            axisDiv = axis;
        }
    }

    double result = (average * numOfPts + meanRecall * numOfRecalls) / (numOfPts + numOfRecalls ) * volume;
    double error = fabs(meanRecall - average ) * volume;
    double tolerance = absAcc + fabs(result ) * relAcc ;

    if (error < tolerance ) {
        return result ;
    }

    double lowerBound_second[dim];
    double upperBound_second[dim];
    for(int axis = 0; axis < dim; axis++) {
        lowerBound_second[axis] = lowerBound[axis] ;
    }

    for(int axis = 0; axis < dim; axis++) {
        upperBound_second[axis] = upperBound[axis];
    }
    lowerBound_second[axisDiv] = (lowerBound[axisDiv] + upperBound[axisDiv]) / 2;
    upperBound_second[axisDiv] = (lowerBound[axisDiv] + upperBound[axisDiv]) / 2 ;

    double resultLeft   = plain_montecarlo_stratifiedSampling(dim, func, lowerBound       , upperBound_second, absAcc / sqrt(2), relAcc, numToTheLeft[axisDiv] , averageLeft[axisDiv]  );
    double resultRight  = plain_montecarlo_stratifiedSampling(dim, func, lowerBound_second, upperBound       , absAcc / sqrt(2), relAcc, numToTheRight[axisDiv], averageRight[axisDiv] );

    return resultLeft + resultRight ;

}
