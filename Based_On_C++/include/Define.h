#ifndef DEFINE_DEFINE
#define DEFINE_DEFINE
#include <iostream>
#include <stdio.h>

#define PI 3.14159265359

typedef double Scalar;

extern Scalar R;
extern Scalar L;
extern Scalar H;
extern int Nx;
extern int Nz;

extern Scalar kx1;
extern Scalar kx2;

void assemble_ls_array(Scalar** ls1, Scalar** ls2, int Nx);

void triplet_2C_standard(int* I, int* J, int nval);


#endif