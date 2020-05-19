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

    IF we hope to restore the original eigenvectors 'x', then we need to SOLVE Q^T*(S*x)*Q = y!
*/


extern "C" {

    extern int dpbstf_(char*, int*, int*, double*, int*, int*);

    extern int dsbgst_(char*, char*, int*, int*, int*, double*, int*, double*, int*, double*, int*, double*, int*);

    extern int dsbtrd_(char*, char*, int*, int*, double*, int*, double*, double*, double*, int*, double*, int*);

    extern int dstevr_(char*, char*, int*, double*, double*, double*, double*, int*, int*, double*, int*, double*, double*, int*, int*, double*, int*, int*, int*, int*);

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
    //std::cout << matrix_A << std::endl;
    printf("Matrix A has bandwidth %d, in Band Storage Format:  \n", kd_a);
    Eigen::MatrixXd MATRIXA = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (A_array, NA_row, NA_col);
    //std::cout << MATRIXA << std::endl;
    printf("Matrix B in Eigen:  \n");
    //std::cout << matrix_B << std::endl;
    printf("Matrix B has bandwidth %d, in Band Storage Format:  \n", kd_b);
    //Eigen::MatrixXd MATRIXB = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (B_array, NB_row, NB_col);
    //std::cout << MATRIXB << std::endl;
    //std::cout << std::endl;
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
    //Eigen::MatrixXd Cholesky_B = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (B_array, NB_row, NB_col);
    //std::cout << Cholesky_B << std::endl; 
    Eigen::SparseMatrix<Scalar> Cholesky_B_SpEigen = Eigen::SparseMatrix<Scalar> (NB_col, NB_col);
    //CholFactorBandToEigenSparse(Cholesky_B_SpEigen, B_array, NB_row, NB_col);
    printf("Cholesky factorization of Matrix B in Eigen format:\n");
    //std::cout << Cholesky_B_SpEigen << std::endl;

    /*
        2. Reduction of GEP:   A*x = lambda*B*x
        Input B should be the banded factor returned by dpbstf, B_array;
    */
   char VECT_IFX = 'V';
   // UPLO: UPLO_STORE
   // N:    N_B
   int KD_A = kd_a;
   // KB:   KD_B
   // AB:   A_array
   // LDAB: NA_row
   // BB:   B_array
   // LDBB: NB_row
   int LDX = N_B;
   double* X = (double*) calloc(LDX*N_B, sizeof(double));   // Here we need to extract the transformation matrix X in order to restore the eigenvectors
   double* WORK_ROTATION2 = (double*) calloc(2*N_B, sizeof(double));
   int INFO_DSBGST = -1;    
   dsbgst_(&VECT_IFX, &UPLO_STORE, &N_B, &KD_A, &KD_B, A_array, &NA_row, B_array, &NB_row, X, &N_B, WORK_ROTATION2, &INFO_DSBGST);
   // print SOLVER information & results
   Eigen::MatrixXd MATRIXX = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (X, LDX, N_B);
   Eigen::SparseMatrix<Scalar> MATRIXC = Eigen::SparseMatrix<Scalar> (N_B, N_B);
   //SymBandPlainCToEigenSparse(MATRIXC, A_array, NA_row, N_B);
   printf("DSBGST SOLVER return INFO:   %d  and Reduced standard form Matrix C in Eigen format:\n", INFO_DSBGST);
   //std::cout << MATRIXC << std::endl;

   char fpath_matrixX[] = "X";
   write_matrix(MATRIXX, fpath_matrixX);

   /*
        3.  Reduces a real symmetric band matrix to symmetric tridiagonal form;
   */
   char VECT_IFQ = 'V'; // Need Q to restore the eigenvectors
   // UPLO: UPLO_STORE
   // N:    N_B
   // KD:   KD_A
   // AB:   A_array, returned by DSBGST subroutine
   // LDAB: NA_row
   double* D = (double*) calloc(N_B, sizeof(double));
   double* E = (double*) calloc(N_B - 1, sizeof(double));
   int LDQ = N_B;
   double* Q = (double*) calloc(LDQ*N_B, sizeof(double));
   double* WORK_ROTATION3 = (double*) calloc(N_B, sizeof(double));
   int INFO_DSBTRD = -1;
   dsbtrd_(&VECT_IFQ, &UPLO_STORE, &N_B, &KD_A, A_array, &NA_row, D, E, Q, &LDQ, WORK_ROTATION3, &INFO_DSBTRD);
   // print SOLVER information & results
   //Eigen::SparseMatrix<Scalar> MATRIXT = Eigen::SparseMatrix<Scalar> (N_B, N_B);
   Eigen::MatrixXd MATRIXQ = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (Q, LDQ, LDQ);
   //TridiagArraysToEigenSparse(MATRIXT, D, E, N_B);
   printf("DSBTRD SOLVER return INFO:   %d  and Reduced tridiagnonal Matrix T in Eigen format:\n", INFO_DSBTRD);
   //std::cout << MATRIXT << std::endl;
   //   Record Matrix Q
   char fpath_matrixQ[] = "Q";
   write_matrix(MATRIXQ, fpath_matrixQ);


   /*
        4.  Compute Eigenvalues and Eigenvectors;
        Some tips for the parameters used in routine DSTEVR:
                a.  If lwork (or liwork) has any of admissible sizes, which is no less than the minimal value described, 
                then the routine completes the task, though probably not so fast as with a recommended workspace, 
                and provides the recommended workspace in the first element of the corresponding array (work, iwork) on exit. 
                Use this value (work(1), iwork(1)) for subsequent runs.
                b.  If lwork = -1 (liwork = -1), then the routine returns immediately and provides the recommended workspace 
                in the first element of the corresponding array (work, iwork). This operation is called a workspace query.
                c.  If lwork (liwork) is less than the minimal required value and is not equal to -1, then the routine returns 
                immediately with an error exit and does not provide any information on the recommended workspace.
   */
   char JOBINFO = 'V';
   char COMPUTEMODE = 'I';
   // N:    N_B
   // D:    D
   // E:    E
   double VL, VU;   // not referenced
   int SMALLEST_EIGVAL_INDEX = 1;
   int LARGEST_EIGVAL_INDEX = 1000;
   double ABSTOL = 1e-10;
   int NUMOFEIG = LARGEST_EIGVAL_INDEX - SMALLEST_EIGVAL_INDEX + 1;
   double* eig_vals = (double*) calloc(N_B, sizeof(double));
   int LDZ = N_B;
   double* eig_vecs = (double*) calloc(LDZ*NUMOFEIG, sizeof(double));
   int* ISUPPZ = (int*) calloc(2*NUMOFEIG, sizeof(int));
   // Temporal storage size, will be re-allocated after workspace querry;
   double* WORK_DSTEVR = (double*) calloc(1, sizeof(double));
   int LWORK_DSTEVR = -1;
   int* IWORK_DSTEVR = (int*) calloc(1, sizeof(int));
   int LIWORK_DSTEVR = -1;
   int INFO_DSTEVR = -1;
   dstevr_(&JOBINFO, &COMPUTEMODE, &N_B, D, E, &VL, &VU, &SMALLEST_EIGVAL_INDEX, &LARGEST_EIGVAL_INDEX, &ABSTOL, &NUMOFEIG, eig_vals, eig_vecs, &LDZ, ISUPPZ, WORK_DSTEVR, &LWORK_DSTEVR, IWORK_DSTEVR, &LIWORK_DSTEVR, &INFO_DSTEVR);
   // print workspace querry information;
   LWORK_DSTEVR = (int)WORK_DSTEVR[0];
   LIWORK_DSTEVR = IWORK_DSTEVR[0];
   printf("Optimal size for LWORK and LIWORK:\t%d\t%d\n", LWORK_DSTEVR, LIWORK_DSTEVR);
   // re-allocate WORK and IWORK array size;
   WORK_DSTEVR = (double*) calloc(LWORK_DSTEVR, sizeof(double));
   IWORK_DSTEVR = (int*) calloc(LIWORK_DSTEVR, sizeof(int));
   dstevr_(&JOBINFO, &COMPUTEMODE, &N_B, D, E, &VL, &VU, &SMALLEST_EIGVAL_INDEX, &LARGEST_EIGVAL_INDEX, &ABSTOL, &NUMOFEIG, eig_vals, eig_vecs, &LDZ, ISUPPZ, WORK_DSTEVR, &LWORK_DSTEVR, IWORK_DSTEVR, &LIWORK_DSTEVR, &INFO_DSTEVR);
   // print SOLVER information & results
   printf("DSTEVR SOLVER return INFO:   %d\n", INFO_DSTEVR);
   Eigen::MatrixXd MATRIX_EIGVECS = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (eig_vecs, LDZ, NUMOFEIG);
   Eigen::MatrixXd MATRIX_EIGVALS = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (eig_vals, NUMOFEIG, 1);
   char fpath_eig_vecs[] = "eig_vecs";
   char fpath_eig_vals[] = "eig_vals";
   write_matrix(MATRIX_EIGVECS, fpath_eig_vecs);
   write_matrix(MATRIX_EIGVALS, fpath_eig_vals);

   /*
        5.  Restore original eigenvectors from Q^T*(S*x)*Q = y
        Q^T*X^(-1)* x = y,  --->    x = X*Q* y
        X is in dsbgst call and Q is in dsbtrd call.
   */
   Eigen::MatrixXd EIGVECS = MATRIXX*MATRIXQ*MATRIX_EIGVECS;


   Eigen::MatrixXd matrix_res(NUMOFEIG, 1);
   for(int i=0; i<NUMOFEIG; i++){
       matrix_res(i) = ((matrix_A - eig_vals[i]*matrix_B)*EIGVECS.col(i)).lpNorm<2>();
   }
   char fpath_residues[] = "eig_res";
   write_matrix(matrix_res, fpath_residues);


}