#include <iostream>
#include "Eigen_Lapack_Interface.h"
#include "DifferenceOperator.h"
#include "FileIO.h"
#include "MatrixAssembly.h"

extern "C" {

extern int dpotrf_(char*,int*,double*,int*,int*);

extern int dsygst_(int*,char*,int*,double*,int*,double*,int*,int*);

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

    /*
        Dense Matrix Test
    */
    Eigen::MatrixXd A(matrix_A);
    Eigen::MatrixXd B(matrix_B);

    std::cout << "Matrix A " << std::endl;
    std::cout << A << std::endl;
    std::cout << "Matrix B " << std::endl;
    std::cout << B << std::endl;
    std::cout << std::endl;


    Scalar* A_array, * B_array;
    int Nr, Nc;
    EigenDenseToPlainC(A, &A_array, &Nr, &Nc);
    EigenDenseToPlainC(B, &B_array, &Nr, &Nc);
    
    std::cout << std::endl;

    char UPLO_CD_B = 'U';
    int INFO_DPOTRF = -666;
    dpotrf_(&UPLO_CD_B, &Nr, B_array, &Nr, &INFO_DPOTRF);
    for(int i=0; i<Nr*Nc; i++){
        printf("%f\t", B_array[i]);
    }

    int EIGENTYPE = 1;
    char UPLO_A = UPLO_CD_B;
    int INFO_DSYGST = -666;
    dsygst_(&EIGENTYPE, &UPLO_A, &Nr, A_array, &Nr, B_array, &Nr, &INFO_DSYGST);

    // print Matrix A and B to have a look
    Eigen::MatrixXd MATRIXA = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (A_array, Nr, Nc);
    Eigen::MatrixXd MATRIXB = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (B_array, Nr, Nc);
    printf("Solver Info:    %d\t%d\n\n", INFO_DPOTRF, INFO_DSYGST);
    
    Eigen::MatrixXd MATRIXBt = MATRIXB.transpose();
    std::cout << MATRIXA << std::endl;
    std::cout << std::endl;
    std::cout << MATRIXB << std::endl;
    std::cout << std::endl;
    //std::cout << MATRIXB*MATRIXBt << std::endl;

}
