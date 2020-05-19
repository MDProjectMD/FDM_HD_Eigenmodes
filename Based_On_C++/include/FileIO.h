#ifndef DEFINE_FILEIO
#define DEFINE_FILEIO
#include "Define.h"
#include <fstream>
#include "Eigen/Dense"
#include "Eigen/Sparse"

void WriteCompactArray(char* fpath, int* I, int* J, Scalar* K, int Nx, int Nz);

// parameters order:	L,H,Nx,Nz
void WriteParameter(char* fpath);

void WriteSlipLength(char* fpath, Scalar* ls1, Scalar* ls2);

void process_monitor(int pid, int N);

void write_matrix(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Mat, char* fpath);

#endif
