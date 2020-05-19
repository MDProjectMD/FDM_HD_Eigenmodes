#include <iostream>
#include "Eigen_Lapack_Interface.h"

#define N 5

// dgeev_ is a symbol in the LAPACK library files
extern "C" {
extern int dpotrf_(char*,int*,double*,int*,int*);
}

/*
Details of Cholesky factorization
   1.77   0.10  -0.51   0.93  -0.41
   0.00   0.88   0.99  -0.84   0.36
   0.00   0.00   1.81  -1.32   0.57
   0.00   0.00   0.00   1.42   0.05
   0.00   0.00   0.00   0.00   1.16
*/

int main(){
    double A[N*N] = {
            3.14,  0.00,  0.00,  0.00,  0.00,
            0.17,  0.79,  0.00,  0.00,  0.00,
           -0.90,  0.83,  4.53,  0.00,  0.00,
            1.65, -0.65, -3.70,  5.32,  0.00,
           -0.72,  0.28,  1.60, -1.37,  1.98
    };
    char UPLO_CHAR = 'U';
    int NA = N;
    int INFO_INT = 0;

    dpotrf_(&UPLO_CHAR, &NA, A, &NA, &INFO_INT);

    // print CD result of matrix A in c storage sequence
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            printf("%f\t", A[j*N + i]);
        }
        printf("\n");
    }
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat(4,4);
    mat   <<  1, 2, 3, 4,
            5, 6, 7, 8,
            9,10,11,12,
            13,14,15,16;
    Scalar* mat_array;
    int Nr, Nc;
    EigenDenseToPlainC(mat, &mat_array, &Nr, &Nc);
    for(int i=0; i<16; i++){
        printf("%f\t", mat_array[i]);
    }
    printf("\n");
}