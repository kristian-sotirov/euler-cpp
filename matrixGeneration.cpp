#include <iostream>
#include <vector>
#include <cmath>
#include "mult.h"
#include "matrixGeneration.h"

std::vector<double> createR(std::vector<double> &b, std::vector<std::vector<double>> &A, std::vector<double> &x) {

	int n = b.size();	
	
	std::vector<double> r = b;
	std::vector<double> Ax = multMatrixVec(A, x);
	
	for(int i = 0; i < n; i++) {
		r[i] -= Ax[i];
	}
	
	return r;
}


std::vector<std::vector<double>> createP(std::vector<std::vector<double>> &M) {
	int size = M.size();
	
	std::vector<std::vector<double>> P = std::vector<std::vector<double>>(size, std::vector<double>(size, 0.0));
	
	for(int i = 0; i < size; i++) {
		if(i != 0) {
			P[i][i-1] = M[i][i-1];
		}
		P[i][i] = M[i][i];
		if (i != size - 1) {
			P[i][i+1] = M[i][i+1];
		}
	}
	
	return P;
}

std::vector<std::vector<double>> createMhK(std::vector<std::vector<double>> &M, double h, std::vector<std::vector<double>> &K) {

	int size = M.size();
	std::vector<std::vector<double>> MhK = M;
	
	for(int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			MhK[i][j] -= h*K[i][j];
		}
	}
	
	return MhK;
}

std::vector<double> createhf(double h, std::vector<double> &f) {
	int size = f.size();
	std::vector<double> hf = f;
	
	for(int i = 0; i < size; i++) {
		hf[i] *= h;
	}
	
	return hf;
}

std::vector<double> createb(std::vector<double> &hf, std::vector<std::vector<double>> &MhK, std::vector<double> &u_k) {

	int size = hf.size();
	std::vector<double> b = multMatrixVec(MhK, u_k);
	
	for(int i = 0; i < size; i++) {
		b[i] += hf[i];
	}
	
	return b;
}



