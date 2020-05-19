#include <iostream>
#include "Eigen_Lapack_Interface.h"
#include "DifferenceOperator.h"
#include "FileIO.h"
#include "MatrixAssembly.h"

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
	Scalar** Cm, ** Lm;
	biharmonic_4th_center(&Cm, DX, DZ);
	laplacian_4th_center(&Lm, DX, DZ);

	int* I_bh, *J_bh, *I_lap, *J_lap;
	Scalar* K_bh, *K_lap;
	
	// 2nd or 4th differs only at the boundary conditions difference scheme.  In the bulk region both take the 4th order difference scheme
	int num_val_bh = assemble_biharmonic_matrix_2nd(Cm, ls1, ls2, &I_bh, &J_bh, &K_bh, Nx, Nz);
	triplet_2C_standard(I_bh, J_bh, num_val_bh); // if using matlab, delete it
	std::vector< Eigen::Triplet<double> > coefficients_A;
	for(int n=0; n<num_val_bh; n++){
		coefficients_A.push_back(Eigen::Triplet<double>(I_bh[n], J_bh[n], K_bh[n]));
	}
	free_biharmonic_matrix(Cm, I_bh, J_bh, K_bh);
	Eigen::SparseMatrix<double> matrix_A = Eigen::SparseMatrix<double> (N, N);	
	matrix_A.setFromTriplets(coefficients_A.begin(), coefficients_A.end());

	int num_val_lap = assemble_laplacian_matrix_2nd(Lm, ls1, ls2, &I_lap, &J_lap, &K_lap, Nx, Nz);
	triplet_2C_standard(I_lap, J_lap, num_val_lap);	// if using matlab, delete it
	std::vector< Eigen::Triplet<double> > coefficients_B;
	for(int n=0; n<num_val_lap; n++){
		coefficients_B.push_back(Eigen::Triplet<double>(I_lap[n], J_lap[n], 0. - K_lap[n]));
	}
	Eigen::SparseMatrix<double> matrix_B = Eigen::SparseMatrix<double> (N, N);
	matrix_B.setFromTriplets(coefficients_B.begin(), coefficients_B.end());
	free_laplacian_matrix(Lm, I_lap, J_lap, K_lap);

    // Transform to Band Matrix A and B in plain C format;
    int kd_b = EigenSymBandWidth(matrix_B);
    int kd_a = EigenSymBandWidth(matrix_A);
    Scalar* B_array, * A_array;             // BAND STORAGE pointers
    int NB_row, NB_col, NA_row, NA_col;     // BAND STORAGE size
    EigenSpSymToBandPlainC(matrix_B, kd_b, &B_array, &NB_row, &NB_col);
    EigenSpSymToBandPlainC(matrix_A, kd_a, &A_array, &NA_row, &NA_col);
    printf("Matrix A has bandwidth %d, Matrix B has bandwidth %d \n", kd_a, kd_b);

    /*********************************  LAPACK SOLVING  Begin   HERE    ***********************************/
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
    char fpath_Cholesky_B[] = "Cholesky_B";
    printf("DPBSTF SOLVER return INFO: %d and prepare writing Band Storage of Matrix B into file %s ......\n", INFO_DPBSTF, fpath_Cholesky_B);
    Eigen::MatrixXd Cholesky_B_Band = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (B_array, NB_row, NB_col);
    write_matrix(Cholesky_B_Band, fpath_Cholesky_B);
    printf("Cholesky factorization of Matrix B in Band Storage format has been written in %s:\n", fpath_Cholesky_B);

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
    Eigen::MatrixXd MATRIXC = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (A_array, NA_row, NA_col);
    Eigen::MatrixXd MATRIXX = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (X, LDX, N_B);
    char fpath_matrixC[] = "Reduced_C";
    printf("DSBGST SOLVER return INFO:   %d  and prepare writing Reduced standard form Matrix C in Band Storage Format into file %s ......\n", INFO_DSBGST, fpath_matrixC);
    write_matrix(MATRIXC, fpath_matrixC);
    printf("Reduced standard form Matrix C in Band Storage Format has been written in %s:\n", fpath_matrixC);
    char fpath_matrixX[] = "X";
    printf("Prepare writing transfomation matrix X into file %s ...... \n", fpath_matrixX);
    write_matrix(MATRIXX, fpath_matrixX);
    printf("Transfomation matrix X has been written in %s:\n", fpath_matrixX);

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
    char fpath_matrixQ[] = "Q";
    char fpath_arrayD[] = "D";
    char fpath_arrayE[] = "E";
    Eigen::MatrixXd MATRIXQ = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (Q, LDQ, LDQ);
    Eigen::MatrixXd ARRAYD = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (D, N_B, 1);
    Eigen::MatrixXd ARRAYE = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (E, N_B - 1, 1);
    printf("DSBTRD SOLVER return INFO:   %d  and prepare writing diagonal D and off-diagonal E arrays into files %s and %s ......:\n", INFO_DSBTRD, fpath_arrayD, fpath_arrayE);
    write_matrix(ARRAYD, fpath_arrayD);
    write_matrix(ARRAYE, fpath_arrayE);
    printf("diagonal D and off-diagonal E arrays has been written into files.\n");
    printf("prepare writing transfomation matrix Q into file %s ......:\n", fpath_matrixQ);
    write_matrix(MATRIXQ, fpath_matrixQ);
    printf("Transfomation matrix Q has been written in %s:\n", fpath_matrixQ);

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
    int LARGEST_EIGVAL_INDEX = 50;
    double ABSTOL = 1e-12;
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
    printf("Optimal size for LWORK and LIWORK are:\t%d\t%d\n", LWORK_DSTEVR, LIWORK_DSTEVR);
    // re-allocate WORK and IWORK array size;
    WORK_DSTEVR = (double*) calloc(LWORK_DSTEVR, sizeof(double));
    IWORK_DSTEVR = (int*) calloc(LIWORK_DSTEVR, sizeof(int));
    dstevr_(&JOBINFO, &COMPUTEMODE, &N_B, D, E, &VL, &VU, &SMALLEST_EIGVAL_INDEX, &LARGEST_EIGVAL_INDEX, &ABSTOL, &NUMOFEIG, eig_vals, eig_vecs, &LDZ, ISUPPZ, WORK_DSTEVR, &LWORK_DSTEVR, IWORK_DSTEVR, &LIWORK_DSTEVR, &INFO_DSTEVR);
    // print SOLVER information & results
    printf("DSTEVR SOLVER return INFO:   %d\n", INFO_DSTEVR);
    char fpath_eig_vecs[] = "eig_vecs";
    char fpath_eig_vals[] = "eigen_values";
    printf("prepare writing eigenvalues and intermediate eigenvectors into file %s and %s ......\n", fpath_eig_vals, fpath_eig_vecs);
    Eigen::MatrixXd MATRIX_EIGVECS = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (eig_vecs, LDZ, NUMOFEIG);
    Eigen::MatrixXd MATRIX_EIGVALS = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (eig_vals, NUMOFEIG, 1);
    write_matrix(MATRIX_EIGVECS, fpath_eig_vecs);
    write_matrix(MATRIX_EIGVALS, fpath_eig_vals);
    printf("eigenvalues and intermediate eigenvectors have been recorded\n");

    /*
        5.  Restore original eigenvectors from Q^T*(S*x)*Q = y
        Q^T*X^(-1)* x = y,  --->    x = X*Q* y
        X is in dsbgst call and Q is in dsbtrd call.
    */
    printf("Begin to restore the original eigenvectors ......\n");
    Eigen::MatrixXd EIGVECS = MATRIXX*MATRIXQ*MATRIX_EIGVECS;
    printf("Begin to check the residue errors ......\n");
    Eigen::MatrixXd matrix_res(NUMOFEIG, 1);
    for(int i=0; i<NUMOFEIG; i++){
        matrix_res(i) = ((matrix_A - eig_vals[i]*matrix_B)*EIGVECS.col(i)).lpNorm<2>();
    }
    char fpath_residues[] = "eig_res";
    write_matrix(matrix_res, fpath_residues);
    char fpath_eig_vecs_restored[] = "eigen_vectors";
    printf("Begin to write the restored eigenvectors into file %s ......", fpath_eig_vecs_restored);
    write_matrix(EIGVECS, fpath_eig_vecs_restored);
    printf("Computation completed!\n");

}