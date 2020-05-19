// self-defined package must be located below the default package !!!
#include "Define.h" 
#include "DifferenceOperator.h"
#include "FileIO.h"
#include "MatrixAssembly.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "Spectra/SymGEigsSolver.h"
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>

int main() {
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
	printf("%f", H);
	printf("%f", L);
	
	// 2nd or 4th differs only at the boundary conditions difference scheme.  In the bulk region both take the 4th order difference scheme
	int num_val_bh = assemble_biharmonic_matrix_2nd(Cm, ls1, ls2, &I_bh, &J_bh, &K_bh, Nx, Nz);
	char fpath_biharmonic[] = "./biharmonic_py.txt";
	triplet_2C_standard(I_bh, J_bh, num_val_bh); // if using matlab, delete it
	//WriteCompactArray(fpath_biharmonic, I_bh, J_bh, K_bh, Nx, Nz);
	std::vector< Eigen::Triplet<double> > coefficients_A;
	for(int n=0; n<num_val_bh; n++){
		coefficients_A.push_back(Eigen::Triplet<double>(I_bh[n], J_bh[n], K_bh[n]));
	}
	free_biharmonic_matrix(Cm, I_bh, J_bh, K_bh);
	// IMPORTANT !!!!!
	// Here only the default ColMajor format is support in this solver;
	// IF RowMajor format is used here, it will report Segmentation fault 11 !!
	Eigen::SparseMatrix<double> matrix_A = Eigen::SparseMatrix<double> (N, N);	
	matrix_A.setFromTriplets(coefficients_A.begin(), coefficients_A.end());

	int num_val_lap = assemble_laplacian_matrix_2nd(Lm, ls1, ls2, &I_lap, &J_lap, &K_lap, Nx, Nz);
	char fpath_laplacian[] = "./laplacian_py.txt";
	triplet_2C_standard(I_lap, J_lap, num_val_lap);	// if using matlab, delete it
	//WriteCompactArray(fpath_laplacian, I_lap, J_lap, K_lap, Nx, Nz);
	std::vector< Eigen::Triplet<double> > coefficients_B;
	for(int n=0; n<num_val_lap; n++){
		coefficients_B.push_back(Eigen::Triplet<double>(I_lap[n], J_lap[n], 0. - K_lap[n]));
	}
	Eigen::SparseMatrix<double> matrix_B = Eigen::SparseMatrix<double> (N, N);
	matrix_B.setFromTriplets(coefficients_B.begin(), coefficients_B.end());
	free_laplacian_matrix(Lm, I_lap, J_lap, K_lap);

	//std::cout << Eigen::MatrixXd(matrix_A) << std::endl;
	//std::cout << Eigen::MatrixXd(matrix_B) << std::endl;

	/* 
	Construct matrix operation object using the wrapper classes
				****** CASE I: 2nd order accuracy for BCs ***********
	In this case, the FDM matrix is symmetric, so adapting the SymGEigsSolver to compute it.
	*/
	Spectra::SparseSymMatProd<double> op(matrix_A);
	Spectra::SparseCholesky<double>  Bop(matrix_B);
	// Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
	unsigned int nev_eig = 2;
	unsigned int ncv_Kry = 6*nev_eig;
	Spectra::SymGEigsSolver<double,  Spectra::SMALLEST_MAGN,  Spectra::SparseSymMatProd<double>,  Spectra::SparseCholesky<double>, Spectra::GEIGS_CHOLESKY> geigs(&op, &Bop, nev_eig, ncv_Kry);	
	// Initialize and compute
	geigs.init();
	int nconv = geigs.compute();
	std::cout << nconv << std::endl;
	Eigen::VectorXd eig_vals;
	Eigen::MatrixXd eig_vecs;
	if(geigs.info() == Spectra::SUCCESSFUL){
		eig_vals = geigs.eigenvalues();
		eig_vecs = geigs.eigenvectors();
		//std::cout << "Eigenvalues: " << eig_vals(0) << std::endl;
		//std::cout << "Error L2 norm: " << ((matrix_A - eig_vals(0)*matrix_B)*eig_vecs.leftCols(1)).lpNorm<2>() << std::endl;
	}else{
		std::cout << "Eigenvectors solver failed!" << std::endl;
	}
	// Error check by compute the norm of (A - \lambda B)*eig_vec
	Eigen::VectorXd eig_res(nev_eig);
	for(int n=0; n<nev_eig; n++){
		eig_res(n) = ((matrix_A - eig_vals(0)*matrix_B)*eig_vecs.col(n)).lpNorm<2>();
	}

	/* 
	Construct matrix operation object using the wrapper classes
				****** CASE II: 4th order accuracy ***********
	In this case, the FDM matrix is asymmetric, so CANNOT adapt the SymGEigsSolver to compute!
	*/



	/*
	char fpath_setting[] = "./setting.txt";
	WriteParameter(fpath_setting);

	char fpath_sliplen[] = "./ls.txt";
	WriteSlipLength(fpath_sliplen, ls1, ls2);	// Nx slip length data points
	*/
	char fpath_eig_vecs[] = "eig_vecs";
	char fpath_eig_vals[] = "eig_vals";
	write_matrix(eig_vecs, fpath_eig_vecs);
	write_matrix(eig_vals, fpath_eig_vals);
}