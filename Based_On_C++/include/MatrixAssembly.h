#ifndef DEFINE_MATASSEMBLY
#define DEFINE_MATASSEMBLY
#include "DifferenceOperator.h"
#include "FileIO.h"

/*
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	Nz-1
	^
	|
	2
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	* * * * * * * * * * * * * * (1 --> Nx-1) * * * * * * * * * * * * * 
	i and j increase along the direction [-L --> L](1~Nx-1) and [-H --> H](2~Nz-1)

	shrink above inner grid points from 2D to 1D sequence according to the j-principal order [(j-2)(Nx-1)+i] 
*/

int assemble_biharmonic_matrix_4th(Scalar** BH, Scalar* ls1, Scalar* ls2, int** I, int** J, Scalar** K, int Nx, int Nz);

int assemble_laplacian_matrix_4th(Scalar** LAP, Scalar* ls1, Scalar* ls2, int** I, int** J, Scalar** K, int Nx, int Nz);

int assemble_biharmonic_matrix_2nd(Scalar** BH, Scalar* ls1, Scalar* ls2, int** I, int** J, Scalar** K, int Nx, int Nz);

int assemble_laplacian_matrix_2nd(Scalar** LAP, Scalar* ls1, Scalar* ls2, int** I, int** J, Scalar** K, int Nx, int Nz);

void free_biharmonic_matrix(Scalar** BH, int* I, int* J, Scalar* K);

void free_laplacian_matrix(Scalar** LAP, int* I, int* J, Scalar* K);

#endif
