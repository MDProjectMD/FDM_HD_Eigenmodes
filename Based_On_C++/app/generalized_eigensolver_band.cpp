#include <iostream>
#include "Eigen_Lapack_Interface.h"
#include "DifferenceOperator.h"
#include "FileIO.h"
#include "MatrixAssembly.h"

/*

    In LAPACK, following routines will be used for computing the Generalized eigenvalue/eigenvector problem:
    dpbstf():   DPBSTF computes a split Cholesky factorization of a real symmetric positive definite band matrix A;
                on exit, the factor S (A = S^T*S) will be returned; Of course in the same format as the intput band matrix A [dimension (LDAB,N)]
                S is a band matrix of the same bandwidth as A.
                http://www.netlib.org/lapack/explore-html/da/dba/group__double_o_t_h_e_rcomputational_ga65461621199058b0eddf6686712c1dd4.html

    dsbgst():   DSBGST reduces a real symmetric-definite banded generalized eigenproblem  A*x = lambda*B*x  to standard form  C*y = lambda*y,
                such that C has the same bandwidth as A;
                B must have been previously factorized as S^T*S by DPBSTF, using a split Cholesky factorization. A is overwritten by C = X^T*A*X, 
                where X = S^(-1)*Q and Q is an orthogonal matrix chosen to preserve the bandwidth of A.
                https://www.netlib.org/lapack/explore-html/da/dba/group__double_o_t_h_e_rcomputational_gae32895eca37df3f73da37145f562e707.html

    dsbtrd():   DSBTRD reduces a real symmetric band matrix A to symmetric tridiagonal form T 
                by an orthogonal similarity transformation  Q^T * A * Q = T.
                http://www.netlib.org/lapack/explore-html/da/dba/group__double_o_t_h_e_rcomputational_ga3e341dd5ee93d0f84ed76aa592552648.html#ga3e341dd5ee93d0f84ed76aa592552648

    dstevr():   DSTEVR computes selected eigenvalues and, optionally, eigenvectors of a real symmetric tridiagonal matrix T;
                http://www.netlib.org/lapack/explore-html/dc/dd2/group__double_o_t_h_e_reigen_ga323734560b8bd052fbc474dc2f0b5605.html

*/


extern "C" {

    extern int dpbstf_(char*, int*, int*, double*, int*, int*);

    extern int dsbgst_(char*, char*, int*, int*, int*, double*, int*, double*, int*, double*, int*, double*, int*);

    extern int dsbtrd_(char*, char*, int*, int*, double*, int*, double*, double*, double*, int*, double*, int*);

}


int main(){
    Scalar DX = 2.*L / (Nx - 1);
	Scalar DZ = 2.*H / (Nz - 1);
	int N = (Nx - 1)*(Nz - 2);

	Scalar* ls1, *ls2;
	assemble_ls_array(&ls1, &ls2, Nx);
	Scalar** C, ** L;
	biharmonic_4th_center(&C, DX, DZ);
	laplacian_4th_center(&L, DX, DZ);

	int* I_bh, *J_bh, *I_lap, *J_lap;
	Scalar* K_bh, *K_lap;
	
	// 2nd or 4th differs only at the boundary conditions difference scheme.  In the bulk region both take the 4th order difference scheme
	int num_val_bh = assemble_biharmonic_matrix_2nd(C, ls1, ls2, &I_bh, &J_bh, &K_bh, Nx, Nz);
	triplet_2C_standard(I_bh, J_bh, num_val_bh); // if using matlab, delete it
	//WriteCompactArray(fpath_biharmonic, I_bh, J_bh, K_bh, Nx, Nz);
	std::vector< Eigen::Triplet<double> > coefficients_A;
	for(int n=0; n<num_val_bh; n++){
		coefficients_A.push_back(Eigen::Triplet<double>(I_bh[n], J_bh[n], K_bh[n]));
	}
	free_biharmonic_matrix(C, I_bh, J_bh, K_bh);
	// IMPORTANT !!!!!
	// Here only the default ColMajor format is support in this solver;
	// IF RowMajor format is used here, it will report Segmentation fault 11 !!
	Eigen::SparseMatrix<double> matrix_A = Eigen::SparseMatrix<double> (N, N);	
	matrix_A.setFromTriplets(coefficients_A.begin(), coefficients_A.end());

	int num_val_lap = assemble_laplacian_matrix_2nd(L, ls1, ls2, &I_lap, &J_lap, &K_lap, Nx, Nz);
	char fpath_laplacian[] = "./laplacian_py.txt";
	triplet_2C_standard(I_lap, J_lap, num_val_lap);	// if using matlab, delete it
	//WriteCompactArray(fpath_laplacian, I_lap, J_lap, K_lap, Nx, Nz);
	std::vector< Eigen::Triplet<double> > coefficients_B;
	for(int n=0; n<num_val_lap; n++){
		coefficients_B.push_back(Eigen::Triplet<double>(I_lap[n], J_lap[n], 0. - K_lap[n]));
	}
	Eigen::SparseMatrix<double> matrix_B = Eigen::SparseMatrix<double> (N, N);
	matrix_B.setFromTriplets(coefficients_B.begin(), coefficients_B.end());
	free_laplacian_matrix(L, I_lap, J_lap, K_lap);

    // Transform to Band Matrix A and B in plain C format;
    int kd_b = EigenSymBandWidth(matrix_B);
    int kd_a = EigenSymBandWidth(matrix_A);

    Scalar* B_array, * A_array;
    int NB_row, NB_col, NA_row, NA_col;     // BAND STORAGE size
    EigenSpSymToBandPlainC(matrix_B, kd_b, &B_array, &NB_row, &NB_col);
    EigenSpSymToBandPlainC(matrix_A, kd_a, &A_array, &NA_row, &NA_col);

    printf("Matrix A in Eigen:  \n");
    std::cout << matrix_A << std::endl;
    printf("Matrix A has bandwidth %d, in Band Storage Format:  \n", kd_a);
    Eigen::MatrixXd MATRIXA = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (A_array, NA_row, NA_col);
    std::cout << MATRIXA << std::endl;
    printf("Matrix B in Eigen:  \n");
    std::cout << matrix_B << std::endl;
    printf("Matrix B has bandwidth %d, in Band Storage Format:  \n", kd_b);
    Eigen::MatrixXd MATRIXB = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (B_array, NB_row, NB_col);
    std::cout << MATRIXB << std::endl;
    std::cout << std::endl;
    /*
        GEP:    A*x = lambda*B*x    
        1. Cholesky factorization of Matrix B
    */
    char UPLO_STORE = 'U';
    int N_B = NB_col;
    int KD_B = kd_b;
    // B_array: pointer of Band Storage of Matrix B in plain C
    int LDAB_B = NB_row;
    int INFO_DPBSTF = -1;
    dpbstf_(&UPLO_STORE, &N_B, &KD_B, B_array, &LDAB_B, &INFO_DPBSTF);
    // print SOLVER information & results
    printf("DPBSTF SOLVER return INFO: %d and Cholesky factorization of Matrix B in Band Storage Format:\n", INFO_DPBSTF);
    Eigen::MatrixXd Cholesky_B = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (B_array, NB_row, NB_col);
    std::cout << Cholesky_B << std::endl; 
    Eigen::SparseMatrix<Scalar> Cholesky_B_SpEigen = Eigen::SparseMatrix<Scalar> (NB_col, NB_col);
    CholFactorBandToEigenSparse(Cholesky_B_SpEigen, B_array, NB_row, NB_col);
    printf("Cholesky factorization of Matrix B in Eigen format:\n");
    std::cout << Cholesky_B_SpEigen << std::endl;

    /*
        2. Reduction of GEP:   A*x = lambda*B*x
        Input B should be the banded factor returned by dpbstf, B_array;
    */
   char VECT_IFX = 'N';
   // UPLO: UPLO_STORE
   // N:    N_B
   int KD_A = kd_a;
   // KB:   KD_B
   // AB:   A_array
   // LDAB: NA_row
   // BB:   B_array
   // LDBB: NB_row
   double* X;   // Here we do not want to extract the transformation matrix X
   int LDX = N_B;   // >= 1 anyway
   double* WORK_ROTATION2 = (double*) calloc(2*N_B, sizeof(double));
   int INFO_DSBGST = -1;    
   dsbgst_(&VECT_IFX, &UPLO_STORE, &N_B, &KD_A, &KD_B, A_array, &NA_row, B_array, &NB_row, X, &N_B, WORK_ROTATION2, &INFO_DSBGST);
   // print SOLVER information & results
   Eigen::SparseMatrix<Scalar> MATRIXC = Eigen::SparseMatrix<Scalar> (N_B, N_B);
   SymBandPlainCToEigenSparse(MATRIXC, A_array, NA_row, N_B);
   printf("DSBGST SOLVER return INFO:   %d  and Reduced standard form Matrix C in Eigen format:\n", INFO_DSBGST);
   std::cout << MATRIXC << std::endl;

   /*
        3.  Reduces a real symmetric band matrix to symmetric tridiagonal form;
   */
   char VECT_IFQ = 'N';
   // UPLO: UPLO_STORE
   // N:    N_B
   // KD:   KD_A
   // AB:   A_array, returned by DSBGST subroutine
   // LDAB: NA_row
   double* D = (double*) calloc(N_B, sizeof(double));
   double* E = (double*) calloc(N_B - 1, sizeof(double));
   double* Q;
   int LDQ = N_B;
   double* WORK_ROTATION3 = (double*) calloc(N_B, sizeof(double));
   int INFO_DSBTRD = -1;
   dsbtrd_(&VECT_IFQ, &UPLO_STORE, &N_B, &KD_A, A_array, &NA_row, D, E, Q, &LDQ, WORK_ROTATION3, &INFO_DSBTRD);
   // print SOLVER information & results
   Eigen::SparseMatrix<Scalar> MATRIXT = Eigen::SparseMatrix<Scalar> (N_B, N_B);
   TridiagArraysToEigenSparse(MATRIXT, D, E, N_B);
   printf("DSBTRD SOLVER return INFO:   %d  and Reduced tridiagnonal Matrix T in Eigen format:\n", INFO_DSBTRD);
   std::cout << MATRIXT << std::endl;


   /*
        3.  Compute Eigenvalues and Eigenvectors
   */
   char JOBINFO = 'V';
   
}