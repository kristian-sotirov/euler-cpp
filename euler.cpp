#include <iostream>
#include <vector>
#include <tuple>
#include "euler.h"
#include <string>
#include <cmath>
#include "mult.h"
#include "matrixGeneration.h"

//THIS FUNCTION IS UNFINISHED!!!!!!!!!!
std::tuple<std::vector<double>, double, std::vector<int>> eulerSolver(std::vector<std::vector<double>> &K, std::vector<std::vector<double>> &M, std::vector<double> &u_0, std::vector<double> &f, int N) {
	
	std::vector<std::vector<double>> invM = invLUfact(M);
	std::vector<std::vector<double>> KinvM = multMatrix(K, invM);
	
	double h = gethmax(KinvM);
	h *= 0.8;
	
	std::vector<std::vector<double>> MhK = createMhK(M, h, K);
	std::vector<double> hf = createhf(h, f);
	
	
	std::vector<double> u_k = u_0;	
	std::vector<double> b = createb(hf, MhK, u_k);
	
	std::vector<std::vector<double>> P = createP(M);
	std::vector<std::vector<double>> Pinv = invLUtridiagfact(P);
	
	std::vector<int> k_values = std::vector<int>(N, 0);
	
	for(int i = 1; i <= N; i++) {
		std::tuple<std::vector<double>, int> step_answer = solveSingleStep(u_k, M, b, Pinv);
		u_k = std::get<0>(step_answer);
		k_values[i] = std::get<1>(step_answer);
		b = createb(hf, MhK, u_k);
	}
	
	std::tuple<std::vector<double>, double, std::vector<int>> return_value = std::make_tuple(u_k, h, k_values);
	return return_value;

}

double gethmax(std::vector<std::vector<double>> &KinvM) {
	double max_eigenvalue = findMaxEigenvalue(KinvM);
	double h_max = 2/max_eigenvalue;
	return h_max;
}

double findMaxEigenvalue(std::vector<std::vector<double>> &KinvM) {

	int n = KinvM.size();
	int num_iter = 10000;
	std::vector<double> eigvec = std::vector<double>(n, 0.0);
	eigvec[0] = 1.0;
	double eigvec_norm = 1.0;
	for(int i = 0; i < num_iter; i++) {
		eigvec = multMatrixVec(KinvM, eigvec);
		eigvec_norm = vecMult(eigvec, eigvec);
		for(int j = 0; j < n; j++) {
			eigvec[j] /= eigvec_norm;
		}
	}
	std::vector<double> eigvec2 = multMatrixVec(KinvM,eigvec);
	double eigval = vecMult(eigvec, eigvec2);
	return eigval * eigval;
}

//THIS FUNCTION IS UNFINISHED!!!!!!!!!
std::tuple<std::vector<double>, int> solveSingleStep(std::vector<double> &u_k, std::vector<std::vector<double>> &A, std::vector<double> &b, std::vector<std::vector<double>> &P) {

	int k = 0;
	std::vector<double> r = createR(b, A, u_k);
	std::vector<double> solution_vector(5, 0.0);
	
	
	
	std::tuple<std::vector<double>, int> step_solution = std::make_tuple(solution_vector, k);
	return step_solution;
}

//THIS FUNCTION IS UNFINISHED!!!!!!!!!
void generateNextStep(std::vector<double> &u_k, double, std::vector<double> &z_k) {

}


std::vector<std::vector<double>> invLUtridiagfact(std::vector<std::vector<double>> P) {
	
	int n = P.size();
	std::vector<std::vector<double>> invP = createIDmatrix(n);
	
	tdLUdecomp(P);
	tdforwardSub(P, invP);
	tdbackwardSub(P, invP);
	
	return invP;
}


void tdLUdecomp(std::vector<std::vector<double>> &P) {
	int n = P.size();
	for (int k = 0; k < n-1; k++) {
		P[k+1][k] /= P[k][k];
		P[k+1][k+1] -= P[k+1][k]*P[k][k+1]; 
	}
}


void tdforwardSub(std::vector<std::vector<double>> L, std::vector<std::vector<double>> &ID) {
	int n = L.size();
	
	for(int columns = 0; columns < n; columns++) {
		for (int i = 1; i < n; i++) {
			ID[i][columns] -= L[i][i-1]*ID[i-1][columns];
		}
	}
}


void tdbackwardSub(std::vector<std::vector<double>> U, std::vector<std::vector<double>> &invP) {
	int n = U.size();
	
	for(int columns = 0; columns < n; columns++) {
		invP[n-1][columns] /= U[n-1][n-1];
		for(int i = n-2; i >= 0; i--) {
			invP[i][columns] -= invP[i+1][columns]*U[i][i+1];
			invP[i][columns] /=(U[i][i]);
		}
	}
}

//This functions find the inverse of a matrix using the LU factorisation.
std::vector<std::vector<double>> invLUfact(std::vector<std::vector<double>> M) {

	int n = M.size();
	std::vector<std::vector<double>> invM = createIDmatrix(n);
	
 	LUdecomp(M);
	forwardSub(M,invM);
 	backwardSub(M, invM);
	
	return invM;
}

std::vector<std::vector<double>> createIDmatrix(int size) {

	std::vector<std::vector<double>> ID = std::vector<std::vector<double>>(size, std::vector<double>(size, 0.0));
	for(int i = 0; i < size; i++) {
		ID[i][i] = 1.0;
	}
	return ID;
}

void LUdecomp(std::vector<std::vector<double>> &M) {

	int n = M.size();
	for (int k = 0; k < n-1; k++) {
		for (int j = k+1; j < n; j++) {
			M[j][k] /= M[k][k];
			for (int i = k+1; i < n; i++) {
				M[j][i] -= M[j][k]*M[k][i];
			}
		}
	}
}

void forwardSub(std::vector<std::vector<double>> L, std::vector<std::vector<double>> &ID) {
	int n = L.size();
	
	for(int IDcolumn = 0; IDcolumn < n; IDcolumn++) {
		for (int i = 1; i < n; i++) {
			for (int j = 0; j < i; j++) {
				ID[i][IDcolumn] -= L[i][j]*ID[j][IDcolumn];
			}
		}
	}
}


void backwardSub(std::vector<std::vector<double>> U, std::vector<std::vector<double>> &invM) {
	
	int n = U.size();
	for(int columns = 0; columns < n; columns++) {
		invM[n-1][columns] /= U[n-1][n-1];
		for(int i = n-2; i >=0; i--) {
			for (int j = n-1; j > i; j--) {
				invM[i][columns] -= U[i][j]*invM[j][columns];
			}
			invM[i][columns] /= U[i][i];
		}
	}
}
