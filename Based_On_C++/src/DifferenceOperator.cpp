#include "DifferenceOperator.h"
#include <cstdlib>

void free_biharmonic_4th_center(Scalar** C) {
	for (int i = 0; i < 5; i++) {
		free(C[i]);
	}
	free(C);
}

void biharmonic_4th_center(Scalar*** C, Scalar DX, Scalar DZ) {
	(*C) = (Scalar**)malloc(5 * sizeof(Scalar*));
	for (int i = 0; i < 5; i++) {
		(*C)[i] = (Scalar*)calloc(5, sizeof(Scalar)); // init 0;
	}
	Scalar DZ2i = 1. / DZ / DZ;
	Scalar DX2i = 1. / DX / DX;
	Scalar DX4i = DX2i * DX2i;
	Scalar DZ4i = DZ2i * DZ2i;
	(*C)[0][2] = DZ4i;
	(*C)[1][1] = 2.*DX2i*DZ2i;	(*C)[1][2] = -4.*DZ4i - 4.*DX2i*DZ2i;	(*C)[1][3] = 2.*DX2i*DZ2i;
	(*C)[2][0] = DX4i;	(*C)[2][1] = -4.*DX4i - 4.*DX2i*DZ2i;	(*C)[2][2] = 6.*DZ4i + 6.*DX4i + 8.*DX2i*DZ2i;	(*C)[2][3] = -4.*DX4i - 4.*DX2i*DZ2i;	(*C)[2][4] = DX4i;
	(*C)[3][1] = 2.*DX2i*DZ2i;	(*C)[3][2] = -4.*DZ4i - 4.*DX2i*DZ2i;	(*C)[3][3] = 2.*DX2i*DZ2i;
	(*C)[4][2] = DZ4i;
	std::cout << "assign 4th order biharmonic difference operator matrix" << std::endl;
}


void free_laplacian_4th_center(Scalar** L) {
	for (int i = 0; i < 5; i++) {
		free(L[i]);
	}
	free(L);
}

void laplacian_4th_center(Scalar*** L, Scalar DX, Scalar DZ) {
	(*L) = (Scalar**)malloc(5 * sizeof(Scalar*));
	for (int i = 0; i < 5; i++) {
		(*L)[i] = (Scalar*)calloc(5, sizeof(Scalar)); // init 0;
	}
	Scalar DZ2i = 1. / DZ / DZ;
	Scalar DX2i = 1. / DX / DX;
	(*L)[0][2] = -1. / 12.*DZ2i;
	(*L)[1][2] = 16. / 12.*DZ2i;
	(*L)[2][0] = -1. / 12.*DX2i;	(*L)[2][1] = 16. / 12.*DX2i;	(*L)[2][2] = -30. / 12.*(DX2i + DZ2i);	(*L)[2][3] = 16. / 12.*DX2i;	(*L)[2][4] = -1. / 12.*DX2i;
	(*L)[3][2] = 16. / 12.*DZ2i;
	(*L)[4][2] = -1. / 12.*DZ2i;
	std::cout << "assign 4th order laplacian difference operator matrix" << std::endl;
}  