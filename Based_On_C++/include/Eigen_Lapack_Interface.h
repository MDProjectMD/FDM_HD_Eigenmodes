#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <iostream>
#include "Define.h"

/*
    Matrix is stored as ColMajor format in fortran;
    Eigen dense matrix object ---> shrinked 1D representation in plain C format;
    Scalar: value type defined in Define.h; 
    (*marray): 1D array, store the elements in Eigen matrix object in ColMajor order; On entry, raw pointer without allocating memory;
    n_row:  return the rows number;
    n_col:  return the columns number;
*/
void EigenDenseToPlainC(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Mat, Scalar** marray, int* n_row, int* n_col);


/*
    Retrieve the real symmetric Eigen Sparse matrix bandwidth;
*/
int EigenSymBandWidth(Eigen::SparseMatrix<Scalar> SpMat);

/*
    Transform real symmetric Eigen Sparse matrix to plain C array in Band Storage Format;
    SpMat:      real symmetric eigen matrix
    kd:         number of subdiagonals or superdiagonals, only the lower or upper triangle needs to be stored; 
                Here only dealing with the upper triangle matrix;
    spmarray:   plain C pointer of band matrix on exit;
    nb_row:     number of rows of band matrix on exit;
    nb_col:     number of columns of band matrix on exit;

    Input A(size N*N), output AB(Band format of A), A is assumed ColMajor here
        AB(kd - (j-i), j) = A(i,j) for  max(0, j - kd) <= i <= j (i - kd =< j <= i + kd)
        Size of AB is (kd+1)*N
        Nonreferenced element in AB are set to be 0 by default;
*/
void EigenSpSymToBandPlainC(Eigen::SparseMatrix<Scalar> SpMat, int kd, Scalar** spmarray, int* nb_row, int* nb_col);


/*
    Restore Sparse Symmetric Matrix in Eigen from Band Storage Format in Plain C;
    bmarray:    plain C symmetric matrix to be restored;
    n_row:      Leading Dimension(LDAB), rows number of matrix bmarray;
    n_col:      N, columns of matrix bmarray;
*/
void SymBandPlainCToEigenSparse(Eigen::SparseMatrix<Scalar>& SpMat, Scalar* bmarray, int n_row, int n_col);


/*
    Retrieve the original matrix in Sparse Eigen format returned by DPBSTF subroutines;
    bmarray:    Split Cholesky factor Matrix returned by DPBSTF in Band Storage Format of Plain C;
    nb_row:     Leading Dimension(LDAB), rows number of matrix bmarray;
    nb_col:     N, columns of matrix bmarray;
*/
void CholFactorBandToEigenSparse(Eigen::SparseMatrix<Scalar>& SpMat, Scalar* bmarray, int nb_row, int nb_col);


/*
    Restore Sparse Symmetric Matrix(Tridiagonal) in Eigen from diag & off-diag arrays;
    D:  diagonal elements of the tridiagonal matrix;
    E:  off-diagonal elements of the tridiagonal matrix;
*/
void TridiagArraysToEigenSparse(Eigen::SparseMatrix<Scalar>& SpMat, Scalar* D, Scalar* E, int N);