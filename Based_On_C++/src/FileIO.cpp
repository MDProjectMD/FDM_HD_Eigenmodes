#include "FileIO.h"
#include <cstdlib>

void WriteCompactArray(char* fpath,int* I, int* J, Scalar* K, int Nx, int Nz) {
	int N = 25 * (Nx - 1)*(Nz - 2);
	FILE* fp = fopen(fpath, "w+");
	if (fp == NULL) {
		std::cout << "wrong file path" << std::endl;
		exit(0);
	}
	else {
		for (int i = 0; i < N; i++) {
			fprintf(fp, "%d\t%d\t%f\n", I[i], J[i], K[i]);
		}
	}
	fclose(fp);
}

void WriteParameter(char* fpath) {
	FILE* fp = fopen(fpath, "w+");
	if (fp != NULL) {
		fprintf(fp, "%f\t%f\t%d\t%d\n", L, H, Nx, Nz);
	}
	fclose(fp);
}

void WriteSlipLength(char* fpath, Scalar* ls1, Scalar* ls2) {
	FILE* fp = fopen(fpath, "w+");
	if (fp != NULL) {
		for (int i = 1; i <= Nx - 1; i++) {
			fprintf(fp, "%f\t%f\n", ls1[i], ls2[i]);
		}
		fprintf(fp, "%f\t%f\n", ls1[1], ls2[1]);
	}
	fclose(fp);
}

void process_monitor(int pid, int N) {
	int percent = (double)pid / N * 100;
	if (percent % 5 == 0) {
		std::cout << "Process " << percent << "%" << std::endl;
	}
}

void write_matrix(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Mat, char* fpath){
	std::ofstream fstream(fpath, std::ios::out | std::ios::trunc);
	fstream << Mat << std::endl;
}