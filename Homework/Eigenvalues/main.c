#include <stdio.h>
#include <time.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "jacobi.h"
#include "extraFuncs.h"

void runtime(int Reps, int startRep, char* my_outputFilename, char* gsl_outputFilename, unsigned int* seed);

int main(int argc, char* argv[]){
    printf("\n\n");
    printf("Part A):\n");
    printf("\n");
    unsigned int seed = time(NULL);
    int dims = 5;

    gsl_matrix* matrix      =  gsl_matrix_alloc(dims, dims);
    gsl_matrix* eigVecMat   =  gsl_matrix_alloc(dims, dims);
    gsl_matrix* eigValMat   =  gsl_matrix_alloc(dims, dims);

    set_data_symmetric(matrix, &seed);
    gsl_matrix_memcpy(eigValMat, matrix);

    print_matrix(dims, matrix, "Before diagonalization:");
    jacobiDiag(eigValMat, eigVecMat);
    
    gsl_matrix* tmp             =  gsl_matrix_alloc(dims, dims);
    gsl_matrix* testIdentity    =  gsl_matrix_alloc(dims, dims);
    gsl_matrix* testDiag        =  gsl_matrix_alloc(dims, dims);
    gsl_matrix* testMat         =  gsl_matrix_alloc(dims, dims);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matrix,    eigVecMat, 0.0, tmp          );
    gsl_blas_dgemm(CblasTrans,   CblasNoTrans, 1.0, eigVecMat, tmp,       0.0, testDiag     );
    gsl_blas_dgemm(CblasTrans,   CblasNoTrans, 1.0, eigVecMat, eigVecMat, 0.0, testIdentity );
    gsl_blas_dgemm(CblasNoTrans, CblasTrans,   1.0, eigValMat, eigVecMat, 0.0, tmp          );
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, eigVecMat, tmp,       0.0, testMat      );


    print_matrix(dims, eigVecMat,    "Eigenvectors V:     ");
    print_matrix(dims, testIdentity, "Check V^T  * V = 1:    ");
    print_matrix(dims, eigValMat,    "Eigenvalues D:      ");
    print_matrix(dims, testDiag,     "Check V^T * A * V = D: ");
    print_matrix(dims, matrix,       "Symmetric matrix A:           ");
    print_matrix(dims, testMat,      "Check V * D * V^T = A: ");


    // Part B
    printf("\n\n");
    printf("Part B):\n");
    printf("\n");

    FILE* myOutputFilestream = fopen(argv[1], "w");

    int divisions = 50;
    printf("Dividing into %i divisions\n", divisions);
    double s = 1.0 / (divisions + 1);
    gsl_matrix* hamiltonian = gsl_matrix_alloc(divisions, divisions);
    for(int id = 0; id < divisions - 1; id++){
        gsl_matrix_set(hamiltonian,   id,     id,     -2 );
        gsl_matrix_set(hamiltonian,   id,   id + 1, 1  );
        gsl_matrix_set(hamiltonian, id + 1, id,     1  );
    }
    gsl_matrix_set(hamiltonian,divisions - 1, divisions - 1, -2);
    gsl_matrix_scale(hamiltonian,-1 / s / s);
    printf("Hamiltonian\n");

    gsl_matrix* eigStates = gsl_matrix_alloc(divisions,divisions);
    jacobiDiag(hamiltonian, eigStates);

    printf("\nEigen-energies:\n");
    printf("#\tCalculated\tExact\n");
    for (int energy = 0; energy < divisions / 3; energy++){
        double exact = M_PI*M_PI*(energy + 1)*(energy + 1);
        double calculated = gsl_matrix_get(hamiltonian, energy, energy);
        printf("%i\t%.5g\t%.5g\n", energy, calculated, exact);
    }

    for(int energy = 0; energy < 3; energy++) {
        fprintf(myOutputFilestream, "0\t0\t");
    }
    fprintf(myOutputFilestream, "\n");
    for(int i = 0; i < divisions; i++){
        fprintf(myOutputFilestream, "%.5g\t", (i + 1.0) / (divisions + 1));
        for(int energy = 0; energy < 3; energy++) {
            fprintf(myOutputFilestream, "%.5g\t", gsl_matrix_get(eigStates, i, energy));
        }
        fprintf(myOutputFilestream, "\n");
    }
    fprintf(myOutputFilestream, "1\t");
    for(int energy = 0; energy < 3; energy++) {
        fprintf(myOutputFilestream, "0\t");
    }
    fprintf(myOutputFilestream, "\n");

    gsl_matrix_free(matrix);
    gsl_matrix_free(eigValMat);
    gsl_matrix_free(eigVecMat);
    gsl_matrix_free(tmp);
    gsl_matrix_free(testDiag);
    gsl_matrix_free(testIdentity);
    gsl_matrix_free(testMat);
    gsl_matrix_free(hamiltonian);

    int Reps = 100;
    int startRep  = 50;
    runtime(Reps, startRep, argv[2], argv[3], &seed);

    return 0;
}

