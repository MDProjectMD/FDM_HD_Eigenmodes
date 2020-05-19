#include "Define.h"
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <math.h>

Scalar Ls1(Scalar x) { // lower boundary
	return 1.2 + sin(kx1*x);
}

Scalar Ls2(Scalar x) { // upper boundary
	return 1.2 + sin(kx2*x);
}

void assemble_ls_array(Scalar** ls1, Scalar** ls2, int Nx) {
	Scalar DX = 2.*L / (Nx - 1);
	(*ls1) = (Scalar*)calloc(Nx, sizeof(Scalar));
	(*ls2) = (Scalar*)calloc(Nx, sizeof(Scalar));
	for (int i = 1; i <= Nx - 1; i++) {
		(*ls1)[i] = Ls1(-L + (i - 1)*DX);
		(*ls2)[i] = Ls2(-L + (i - 1)*DX);
	}
}


Scalar R = 0.41;
Scalar L = 0.5;
Scalar H = 0.5;
int Nx = 7;
int Nz = 7;
Scalar kx1 = 0;//PI / L;
Scalar kx2 = 0;//PI / L;

void triplet_2C_standard(int* I, int* J, int nval) {
	// both row and column index should be subtracted by 1
	for (int n = 0; n < nval; n++) {
		I[n] = I[n] - 1;
		J[n] = J[n] - 1;
	}
}

