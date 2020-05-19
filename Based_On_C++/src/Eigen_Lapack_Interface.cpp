#include "Eigen_Lapack_Interface.h"

void EigenDenseToPlainC(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Mat, Scalar** marray, int* n_row, int* n_col){
    Scalar* ptrtmp = NULL;
    if(Mat.IsRowMajor){
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatCol= Mat;
        ptrtmp = MatCol.data();
    }else{
        ptrtmp = Mat.data();
    }

    (*n_row) = Mat.rows();
    (*n_col) = Mat.cols();
    int N = (*n_row)*(*n_col);
    (*marray) = (Scalar*) calloc(N, sizeof(Scalar));
    // memcpy elements to marray;
    if((*marray) != NULL){
        memcpy((*marray), ptrtmp, sizeof(Scalar)*N);
    }else{
        printf("Error allocating memory!\n");   exit(0);
    }
}


int EigenSymBandWidth(Eigen::SparseMatrix<Scalar> SpMat){
    int N_Spmat = -1;
    if(SpMat.rows() == SpMat.cols()){
        N_Spmat = SpMat.rows();
    }else{
        printf("Error calling 'EigenSymBandToPlainC' function for non square matrix! \n");  exit(0);
    }
    int bandwidth = 0;
    for(int k=0; k<SpMat.outerSize(); k++){
        for(Eigen::SparseMatrix<Scalar>::InnerIterator it(SpMat,k); it; ++it){
            int rowIdx = it.row();
            int colIdx = it.col();
            int diff = colIdx - rowIdx;
            if(diff > bandwidth && it.value() != 0){ // To avoid some cases where 
                bandwidth = diff;
            }
        }
    }
    return bandwidth;
}


void EigenSpSymToBandPlainC(Eigen::SparseMatrix<Scalar> SpMat, int kd, Scalar** spmarray, int* n_row, int* n_col){
    int N_Spmat = -1;
    if(SpMat.rows() == SpMat.cols()){
        N_Spmat = SpMat.rows();
    }else{
        printf("Error calling 'EigenSymBandToPlainC' function for non square matrix! \n");  exit(0);
    }
    (*n_row) = kd + 1;
    (*n_col) = N_Spmat;
    int N = (*n_row)*(*n_col);
    (*spmarray) = (Scalar*) calloc(N, sizeof(Scalar));
    for(int j=0; j<N_Spmat; j++){
        int imin = j - kd;
        if(imin <= 0){
            imin = 0;
        }
        for(int i=imin; i<=j; i++){
            int rowIdx = kd - (j-i);
            int colIdx = j;
            (*spmarray)[colIdx*(*n_row) + rowIdx] = SpMat.coeff(i, j);
        }
    }
}

void CholFactorBandToEigenSparse(Eigen::SparseMatrix<Scalar>& SpMat, Scalar* bmarray, int nb_row, int nb_col){
    unsigned int kd = nb_row - 1;
    int m = (nb_col + kd)/(int)2;   // in LAPACK source code, it uses this to find the splitting point m;
    std::vector< Eigen::Triplet<double> > coefficients;
    // retrieve the upper triangle matrix;
    for(int j=0; j<m; j++){
        int imin = j - kd;
        if(imin < 0){
            imin = 0;
        }
        for(int i=imin; i<=j; i++){
            int rowIdxAB = kd + i - j;
            int colIdxAB = j;
            Scalar value = bmarray[colIdxAB*nb_row + rowIdxAB];
            coefficients.push_back(Eigen::Triplet<double>(i, j, value));
        }
    }
    // retrieve the remaining matrix elements;
    for(int j=m; j<nb_col; j++){
        for(int i=0; i<=kd; i++){
            int rowIdxAB = kd - i;
            int colIdxAB = j;
            Scalar value = bmarray[colIdxAB*nb_row + rowIdxAB];
            coefficients.push_back(Eigen::Triplet<double>(colIdxAB, colIdxAB - i, value));
        }
    }
    SpMat.setFromTriplets(coefficients.begin(), coefficients.end());
}

void SymBandPlainCToEigenSparse(Eigen::SparseMatrix<Scalar>& SpMat, Scalar* bmarray, int n_row, int n_col){
    std::vector< Eigen::Triplet<double> > coefficients;
    unsigned int kd = n_row - 1;
    for(int j=0; j<n_col; j++){
        int imin = j - kd;
        if(imin < 0){
            imin = 0;
        }
        for(int i=imin; i<=j; i++){
            int rowIdx = kd + i - j;
            int colIdx = j;
            Scalar value = bmarray[colIdx*n_row + rowIdx];
            coefficients.push_back(Eigen::Triplet<Scalar>(i, j, value));
            if(i!=j){
                coefficients.push_back(Eigen::Triplet<Scalar>(j, i, value));
            }
        }
    }
    SpMat.setFromTriplets(coefficients.begin(), coefficients.end());
}

void TridiagArraysToEigenSparse(Eigen::SparseMatrix<Scalar>& SpMat, Scalar* D, Scalar* E, int N){
    std::vector< Eigen::Triplet<double> > coefficients;
    for(int i=0; i<N-1; i++){
        coefficients.push_back(Eigen::Triplet<Scalar>(i, i+1, E[i]));
    }
    for(int i=1; i<N; i++){
        coefficients.push_back(Eigen::Triplet<Scalar>(i, i-1, E[i-1]));
    }
    for(int i=0; i<N; i++){
        coefficients.push_back(Eigen::Triplet<Scalar>(i, i, D[i]));
    }
    SpMat.setFromTriplets(coefficients.begin(), coefficients.end());
}