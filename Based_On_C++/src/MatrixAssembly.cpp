#include "MatrixAssembly.h"
#include <cstdlib>

void allocate_biharmonic_matrix(int** I, int** J, Scalar** K, int Nx, int Nz) {
	// maximum 25 elements per Phi point
	int N = (Nx - 1)*(Nz - 2);
	(*I) = (int*)calloc(25 * N, sizeof(int));
	(*J) = (int*)calloc(25 * N, sizeof(int));
	(*K) = (Scalar*)calloc(25 * N, sizeof(Scalar));
}
            

int Shift(int i, int Nx) {
	if (i < 1) {
		return i + Nx - 1;
	}
	else if (i > Nx - 1) {
		return i - (Nx - 1);
	}
	else {
		return i;
	}
}

// i = [-2,-1,0,1,2] and j = [-2,-1,0,1,2] refer to the offset-index
Scalar mat_C(Scalar** C, int i, int j) {
	return C[2 + i][2 + j];
}

int shrinkTo1D(int i, int j, int Nx) {
	return (j - 2)*(Nx - 1) + i;
}

void free_biharmonic_matrix(Scalar** BH, int* I, int* J, Scalar* K) {
	free(I);
	free(J);
	free(K);
	for (int i = 0; i < 5; i++) {
		free(BH[i]);
	}
	free(BH);
}

int assemble_biharmonic_matrix_4th(Scalar** BH, Scalar* ls1, Scalar* ls2, int** I, int** J, Scalar** K, int Nx, int Nz) {
	std::cout << "assemble biharmonic matrix in 4th order" << std::endl;
	Scalar DX = 2.*L / (Nx - 1);
	Scalar DZ = 2.*H / (Nz - 1);
	allocate_biharmonic_matrix(I, J, K, Nx, Nz);
	// scan the inner grid points [1~Nx-1]*[2~Nz-1]
	int gridIdx = 0; // record the non-zero C element index in vector I,J and K 
	for (int i = 1; i <= Nx - 1; i++) { // j = 2
		int j = 2;
		int I_idx = i;
		double R1 = (18.*DZ - 6.*ls1[i]) / (3.*DZ + 11.*ls1[i]);
		double R2 = -(4.*ls1[i] + 6.*DZ) / (11.*ls1[i] + 3.*DZ);
		double R3 = (ls1[i] + DZ) / (11.*ls1[i] + 3.*DZ);
		// loop for neighbor grid points
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = 0; beta <= 2; beta++) {
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = j + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(BH, alpha, beta);
				if (alpha == 0) {
					if (beta == 0) {
						(*K)[gridIdx] += R1 * mat_C(BH, 0, -2);
					}
					else if (beta == 1) {
						(*K)[gridIdx] += R2 * mat_C(BH, 0, -2);
					}
					else if (beta == 2) { 
						(*K)[gridIdx] += R3 * mat_C(BH, 0, -2);
					}
				}
				gridIdx++;
			}
		}
	}

	for (int i = 1; i <= Nx - 1; i++) { // j = 3
		int j = 3;
		int I_idx = Nx - 1 + i;
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = -1; beta <= 2; beta++) {
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = j + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(BH, alpha, beta);
				gridIdx++;
			}
		}
	}
	// grid points in bulk region
	for (int j = 4; j <= Nz - 3; j++) {
		process_monitor(j, Nz);
		for (int i = 1; i <= Nx - 1; i++) {
			int I_idx = shrinkTo1D(i, j, Nx);
			for (int alpha = -2; alpha <= 2; alpha++) {
				for (int beta = -2; beta <= 2; beta++) {
					int i_shift = Shift(i + alpha, Nx);
					int j_shift = j + beta;
					int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
					(*I)[gridIdx] = I_idx;
					(*J)[gridIdx] = J_idx;
					(*K)[gridIdx] = mat_C(BH, alpha, beta);
					gridIdx++;
				}
			}
		}
	}

	for (int i = 1; i <= Nx - 1; i++) {
		int j = Nz - 2;
		int I_idx = shrinkTo1D(i, j, Nx);
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = -2; beta <= 1; beta++) {
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = j + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(BH, alpha, beta);
				gridIdx++;
			}
		}
	}

	for (int i = 1; i <= Nx - 1; i++) {
		int j = Nz - 1;
		int I_idx = shrinkTo1D(i, j, Nx);
		double R1 = (18.*DZ - 6.*ls2[i]) / (3.*DZ + 11.*ls2[i]);
		double R2 = -(4.*ls2[i] + 6.*DZ) / (11.*ls2[i] + 3.*DZ);
		double R3 = (ls2[i] + DZ) / (11.*ls2[i] + 3.*DZ);
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = -2; beta <= 0; beta++) {
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = j + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(BH, alpha, beta);
				if (alpha == 0) {
					if (beta == -2) {
						(*K)[gridIdx] += R3 * mat_C(BH, 0, 2);
					}
					else if (beta == -1) {
						(*K)[gridIdx] += R2 * mat_C(BH, 0, 2);
					}
					else if (beta == 0) {
						(*K)[gridIdx] += R1 * mat_C(BH, 0, 2);
					}
				}
				gridIdx++;
			}
		}
	}
	return gridIdx;
}

int assemble_biharmonic_matrix_2nd(Scalar** BH, Scalar* ls1, Scalar* ls2, int** I, int** J, Scalar** K, int Nx, int Nz) {
	std::cout << "assemble biharmonic matrix in 2nd order" << std::endl;
	Scalar DX = 2.*L / (Nx - 1);
	Scalar DZ = 2.*H / (Nz - 1);
	allocate_biharmonic_matrix(I, J, K, Nx, Nz);
	// scan the inner grid points [1~Nx-1]*[2~Nz-1]
	int gridIdx = 0; // record the non-zero C element index in vector I,J and K 
	for (int i = 1; i <= Nx - 1; i++) { // j = 2
		int j = 2;
		int I_idx = i;
		double R = (DZ - 2.*ls1[i]) / (DZ + 2.*ls1[i]);
		// loop for neighbor grid points
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = 0; beta <= 2; beta++) {
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = j + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(BH, alpha, beta);
				if (alpha == 0) {
					if (beta == 0) {
						(*K)[gridIdx] += R * mat_C(BH, 0, -2);
					}
				}
				gridIdx++;
			}
		}
	}

	for (int i = 1; i <= Nx - 1; i++) { // j = 3
		int j = 3;
		int I_idx = Nx - 1 + i;
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = -1; beta <= 2; beta++) {
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = j + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(BH, alpha, beta);
				gridIdx++;
			}
		}
	}
	// grid points in bulk region
	for (int j = 4; j <= Nz - 3; j++) {
		process_monitor(j, Nz);
		for (int i = 1; i <= Nx - 1; i++) {
			int I_idx = shrinkTo1D(i, j, Nx);
			for (int alpha = -2; alpha <= 2; alpha++) {
				for (int beta = -2; beta <= 2; beta++) {
					int i_shift = Shift(i + alpha, Nx);
					int j_shift = j + beta;
					int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
					(*I)[gridIdx] = I_idx;
					(*J)[gridIdx] = J_idx;
					(*K)[gridIdx] = mat_C(BH, alpha, beta);
					gridIdx++;
				}
			}
		}
	}

	for (int i = 1; i <= Nx - 1; i++) {
		int j = Nz - 2;
		int I_idx = shrinkTo1D(i, j, Nx);
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = -2; beta <= 1; beta++) {
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = j + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(BH, alpha, beta);
				gridIdx++;
			}
		}
	}

	for (int i = 1; i <= Nx - 1; i++) {
		int j = Nz - 1;
		int I_idx = shrinkTo1D(i, j, Nx);
		double R = (DZ - 2.*ls2[i]) / (DZ + 2.*ls2[i]);
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = -2; beta <= 0; beta++) {
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = j + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(BH, alpha, beta);
				if (alpha == 0) {
					if (beta == 0) {
						(*K)[gridIdx] += R * mat_C(BH, 0, 2);
					}
				}
				gridIdx++;
			}
		}
	}
	return gridIdx;
}

void allocate_laplacian_matrix(int** I, int** J, Scalar** K, int Nx, int Nz) {
	int N = (Nx - 1)*(Nz - 2);
	(*I) = (int*)calloc(25 * N, sizeof(int));
	(*J) = (int*)calloc(25 * N, sizeof(int));
	(*K) = (Scalar*)calloc(25 * N, sizeof(Scalar));
}

void free_laplacian_matrix(Scalar** LAP, int* I, int* J, Scalar* K) {
	free(I);
	free(J);
	free(K);
	for (int i = 0; i < 5; i++) {
		free(LAP[i]);
	}
	free(LAP);
}

int assemble_laplacian_matrix_2nd(Scalar** LAP, Scalar* ls1, Scalar* ls2, int** I, int** J, Scalar** K, int Nx, int Nz) {
	std::cout << "assemble laplacian matrix in 2nd order" << std::endl;
	Scalar DX = 2.*L / (Nx - 1);
	Scalar DZ = 2.*H / (Nz - 1);
	int gridIdx = 0;
	allocate_laplacian_matrix(I, J, K, Nx, Nz);
	for (int i = 1; i <= Nx - 1; i++) { // j = 2
		int I_idx = shrinkTo1D(i, 2, Nx);
		double R = (DZ - 2.*ls1[i]) / (DZ + 2.*ls1[i]);
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = 0; beta <= 2; beta++) {
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = 2 + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(LAP, alpha, beta);
				if (alpha == 0 && beta == 0) {
					(*K)[gridIdx] += R * mat_C(LAP, 0, -2);
				}
				gridIdx++;
			}
		}

	}

	for (int i = 1; i <= Nx - 1; i++) { // j = 3
		int I_idx = shrinkTo1D(i, 3, Nx);
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = -1; beta <= 2; beta++) {
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = 3 + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(LAP, alpha, beta);
				gridIdx++;
			}
		}
	}

	for (int j = 4; j <= Nz - 3; j++) {
		process_monitor(j, Nz);
		for (int i = 1; i <= Nx - 1; i++) {
			int I_idx = shrinkTo1D(i, j, Nx);
			for (int alpha = -2; alpha <= 2; alpha++) {
				for (int beta = -2; beta <= 2; beta++) {
					int i_shift = Shift(i + alpha, Nx);
					int j_shift = j + beta;
					int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
					(*I)[gridIdx] = I_idx;
					(*J)[gridIdx] = J_idx;
					(*K)[gridIdx] = mat_C(LAP, alpha, beta);
					gridIdx++;
				}
			}
		}
	}

	for (int i = 1; i <= Nx - 1; i++) {
		int j = Nz - 2;
		int I_idx = shrinkTo1D(i, j, Nx);
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = -2; beta <= 1; beta++) { 
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = j + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(LAP, alpha, beta);
				gridIdx++;
			}
		}
	}

	for (int i = 1; i <= Nx - 1; i++) {
		int j = Nz - 1;
		int I_idx = shrinkTo1D(i, j, Nx);
		double R = (DZ - 2.*ls2[i]) / (DZ + 2.*ls2[i]);
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = -2; beta <= 0; beta++) {
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = j + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(LAP, alpha, beta);
				if (alpha == 0) {
					if (beta == 0) {
						(*K)[gridIdx] += R * mat_C(LAP, 0, 2);
					}
				}
				gridIdx++;
			}
		}
	}
	return gridIdx;
}

int assemble_laplacian_matrix_4th(Scalar** LAP, Scalar* ls1, Scalar* ls2, int** I, int** J, Scalar** K, int Nx, int Nz) {
	std::cout << "assemble laplacian matrix in 4th order" << std::endl;
	Scalar DX = 2.*L / (Nx - 1);
	Scalar DZ = 2.*H / (Nz - 1);
	int gridIdx = 0;
	allocate_laplacian_matrix(I, J, K, Nx, Nz);
	for (int i = 1; i <= Nx - 1; i++) { // j = 2
		double R1 = (18.*DZ - 6.*ls1[i]) / (3.*DZ + 11.*ls1[i]);
		double R2 = -(4.*ls1[i] + 6.*DZ) / (11.*ls1[i] + 3.*DZ);
		double R3 = (ls1[i] + DZ) / (11.*ls1[i] + 3.*DZ);
		int I_idx = shrinkTo1D(i, 2, Nx);
		double R = (DZ - 2.*ls1[i]) / (DZ + 2.*ls1[i]);
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = 0; beta <= 2; beta++) {
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = 2 + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(LAP, alpha, beta);
				if (alpha == 0) {
					if (beta == 0) {
						(*K)[gridIdx] += R1 * mat_C(LAP, 0, -2);
					}
					else if (beta == 1) {
						(*K)[gridIdx] += R2 * mat_C(LAP, 0, -2);
					}
					else if (beta == 2) {
						(*K)[gridIdx] += R3 * mat_C(LAP, 0, -2);
					}
					
				}
				gridIdx++;
			}
		}
	}

	for (int i = 1; i <= Nx - 1; i++) { // j = 3
		int I_idx = shrinkTo1D(i, 3, Nx);
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = -1; beta <= 2; beta++) {
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = 3 + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(LAP, alpha, beta);
				gridIdx++;
			}
		}
	}

	for (int j = 4; j <= Nz - 3; j++) {
		process_monitor(j, Nz);
		for (int i = 1; i <= Nx - 1; i++) {
			int I_idx = shrinkTo1D(i, j, Nx);
			for (int alpha = -2; alpha <= 2; alpha++) {
				for (int beta = -2; beta <= 2; beta++) {
					int i_shift = Shift(i + alpha, Nx);
					int j_shift = j + beta;
					int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
					(*I)[gridIdx] = I_idx;
					(*J)[gridIdx] = J_idx;
					(*K)[gridIdx] = mat_C(LAP, alpha, beta);
					gridIdx++;
				}
			}
		}
	}

	for (int i = 1; i <= Nx - 1; i++) {
		int j = Nz - 2;
		int I_idx = shrinkTo1D(i, j, Nx);
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = -2; beta <= 1; beta++) {
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = j + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(LAP, alpha, beta);
				gridIdx++;
			}
		}
	}

	for (int i = 1; i <= Nx - 1; i++) {
		int j = Nz - 1;
		int I_idx = shrinkTo1D(i, j, Nx);
		double R1 = (18.*DZ - 6.*ls2[i]) / (3.*DZ + 11.*ls2[i]);
		double R2 = -(4.*ls2[i] + 6.*DZ) / (11.*ls2[i] + 3.*DZ);
		double R3 = (ls2[i] + DZ) / (11.*ls2[i] + 3.*DZ);
		for (int alpha = -2; alpha <= 2; alpha++) {
			for (int beta = -2; beta <= 0; beta++) {
				int i_shift = Shift(i + alpha, Nx);
				int j_shift = j + beta;
				int J_idx = shrinkTo1D(i_shift, j_shift, Nx);
				(*I)[gridIdx] = I_idx;
				(*J)[gridIdx] = J_idx;
				(*K)[gridIdx] = mat_C(LAP, alpha, beta);
				if (alpha == 0) {
					if (beta == 0) {
						(*K)[gridIdx] += R1 * mat_C(LAP, 0, 2);
					}
					else if (beta == -1) {
						(*K)[gridIdx] += R2 * mat_C(LAP, 0, 2);
					}
					else if (beta == -2) {
						(*K)[gridIdx] += R3 * mat_C(LAP, 0, 2);
					}
				}
				gridIdx++;
			}
		}
	}
	return gridIdx;
}