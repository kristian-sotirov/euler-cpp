#include <iostream>
#include <vector>
#include <tuple>
#include "euler.h"
#include <string>
#include<cmath>

//THIS FUNCTION IS UNFINISHED!!!!!!!!!!
std::tuple<std::vector<double>, double, std::vector<double>> eulerSolver(std::vector<std::vector<double>> &K, std::vector<std::vector<double>> &M, std::vector<double> &u_0, std::vector<double> &f) {
	
	std::tuple<std::vector<double>, double, std::vector<double>> return_value = std::make_tuple(u_0, 0.0, f);
	
	std::vector<std::vector<double>> invM = invLUfact(M);
	std::vector<std::vector<double>> KinvM = multMatrix(K, M);
	
	double h = gethmax(KinvM);
	h *= 0.8;
	
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

double vecMult(std::vector<double> &u, std::vector<double> &v) {
	
	int p = u.size();
	double  mult_value = 0.0;
	
	if (p != v.size()) {
		std::cout << "The vectors cannot be multiplied" << std::endl;
		return mult_value;
	}
	
	for (int i = 0; i < p; i++) {
		mult_value += u[i] * v[i];
	}
	return sqrt(mult_value); 
}

//THIS FUNCTION IS UNFINISHED!!!!!!!!!
std::vector<double> solveSingleStep(std::vector<double> &u_k, std::vector<std::vector<double>> &P) {
	std::vector<double> solution_vector(5, 0.0);
	return solution_vector;
}
//THIS FUNCTION IS UNFINISHED!!!!!!!!!
std::vector<std::vector<double>> createP(std::vector<std::vector<double>> &M) {
	std::vector<std::vector<double>> P = {{0.0,0.0,0.0}, {0.0,0.0,0.0}};
	return P;
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

	std::vector<std::vector<double>> M_decompLU = M;
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


std::vector<std::vector<double>> multMatrix(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B) {

	int m = A.size();
	int p = A[0].size();
	
	if(p != B.size()) {
		std::cout << "The matrices cannot be multiplied" << std::endl;
		return A;
	}
	int n = B[0].size(); 
	
	std::vector<std::vector<double>> C(m, std::vector<double>(n));
	
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < n; j++) {
		
			for(int k = 0; k < p; k++) {
				
				C[i][j] += A[i][k]*B[k][j];
			}
		}
	}
	
	return C;
}

std::vector<double> multMatrixVec(std::vector<std::vector<double>> &A, std::vector<double> &v) {
	int m = A.size();
	int p = A[0].size();
	
	if(p != v.size()) {
		std::cout << "The matrix and the vector cannot be multiplied" << std::endl;
		return v;
	}
	
	std::vector<double> u(m);
	
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < p; j++) {
			u[i] += A[i][j]*v[j];
		}
	}
	return u;
}

std::vector<double> multVecMatrix(std::vector<double> &v, std::vector<std::vector<double>> &B) {
	int p = v.size();
	
	if(p != B.size()) {
		std::cout << "The matrix and the vector cannot be multiplied" << std::endl;
		return v;
	}
	
	int n = B[0].size();
	
	std::vector<double> u(n);
	
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < p; j++) {
			u[i] += v[j]*B[j][i];
		}
	}
	return u;
}
