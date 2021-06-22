#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "jacobi.h"


void jacobiMultiply_right(gsl_matrix* matrix, int firstId, int secondId, double angle){
    /*
     * Multiplies the jacobi matrix JacobiMatrix(firstId, secondId, angle) ( J(p, q, angle) )
     * on the real symmetric matrix matrix (A) from the right, and sets matrix equal to that
     * product; A <-- A * J(p, q, angle)
     *
     *      ¤ matrix    : gsl_matrix* to a real symmetric matrix, A
     *      ¤ firstId   : the first index to rotate in, p
     *      ¤ secondId  : the second index to rotate in, q
     *      ¤ angle     : the angle to rotate by, angle
     *
     */

    double c = cos(angle);
    double s = sin(angle);
    for(int row = 0; row < (matrix -> size1 ); row++){
        double new_aip = c * gsl_matrix_get(matrix, row, firstId)  -  s * gsl_matrix_get(matrix, row, secondId);
        double new_aiq = s * gsl_matrix_get(matrix, row, firstId)  +  c * gsl_matrix_get(matrix, row, secondId);

        gsl_matrix_set(matrix, row, firstId,  new_aip);
        gsl_matrix_set(matrix, row, secondId, new_aiq);
    }
}


void jacobiMultiply_left(gsl_matrix* matrix, int firstId, int secondId, double angle){
    /*
     * Multiplies the jacobi matrix JacobiMatrix(firstId, secondId, angle) ( J(p, q, angle) )
     * on the real symmetric matrix matrix (A) from the left, and sets matrix equal to that
     * product; A <-- J(p, q, angle) * A
     *
     *      ¤ matrix    : gsl_matrix* to a real symmetric matrix, A
     *      ¤ firstId   : the first index to rotate in, p
     *      ¤ secondId  : the second index to rotate in, q
     *      ¤ angle     : the angle to rotate by, angle
     *
     */

    double c = cos(angle);
    double s = sin(angle);
    for(int col = 0; col < (matrix -> size2) ; col++){
        double new_apj =   c * gsl_matrix_get(matrix, firstId, col)  +  s * gsl_matrix_get(matrix, secondId, col);
        double new_aqj = - s * gsl_matrix_get(matrix, firstId, col)  +  c * gsl_matrix_get(matrix, secondId, col);
        gsl_matrix_set(matrix, firstId,  col, new_apj);
        gsl_matrix_set(matrix, secondId, col, new_aqj);
    }
}


void jacobiDiag (gsl_matrix* matrix, gsl_matrix* eigVecMat) {
    /*
     *  Computes the eigenvalue decomposition, using the jacobi diagonalization algorithm
     *  in which elements are jacobi rotated with cyclic sweeps to diagonalize the matrix
     *  (A), which is a real and symmetric matrix. The eigenvalue decomposition consits of
     *
     *      A = V^T * D * V
     *
     *  Where V is an orthonormal matrix consisting of columns that are the eigenvectors of
     *  the matrix A, and D is a diagonal matrix with eigenvalues of A on the diagonal.
     *
     *          ¤ matrix : gsl_matrix* to a real and symmetric matrix, A
     *          ¤ eigVecMat : gsl_matrix* to a square matrix that will hold the matrix
     *                        of eigenvectors, V
     *          ¤ eigValMat : gsl_matrix* to a square matrix that will hold the matrix
     *                        of eigenvalues, D
     *
     */

    // Start by setting V equal to identity matrix in order to do
    // V = I * J_1 * J_2 * .... * J_n, for n sweeps.
    gsl_matrix_set_identity(eigVecMat);

    int dims = matrix -> size1;

    int changed;
    do {
        changed = 0;
        for (int firstId = 0; firstId < dims - 1; firstId++){
            for (int secondId = firstId + 1; secondId < dims; secondId++) {
                double apq = gsl_matrix_get(matrix, firstId,  secondId);
                double app = gsl_matrix_get(matrix, firstId,  firstId );
                double aqq = gsl_matrix_get(matrix, secondId, secondId);
                double angle = 0.5 * atan2(2 * apq, aqq - app);

                double c = cos(angle);
                double s = sin(angle);
                double new_app = c * c * app - 2 * s * c * apq + s * s * aqq;
                double new_aqq = s * s * app + 2 * s * c * apq + c * c * aqq;
                if (new_app != app || new_aqq != aqq){ // do rotation
                    changed = 1;
                    jacobiMultiply_right(matrix,    firstId, secondId,  angle);
                    jacobiMultiply_left( matrix,    firstId, secondId, -angle); // A <- J^T * A * J
                    jacobiMultiply_right(eigVecMat, firstId, secondId,  angle); // V <- V * J
                }
            }
        }
    } while (changed != 0);
}
