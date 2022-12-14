#include <iostream>
#include <vector>
#include <tuple>
#include "euler.h"
#include <string>
#include <cmath>
#include "mult.h"
#include "matrixGeneration.h"

//The function takes as an argument the functions in the following order:
//1) The matrix "K", representing the linear transformation of U.
//2) The matrix "M" representing the discretisation of the derivative.
//3) The vector "U_0", representing the initial condition.
//4) The vector "f", representing the non-linear factor.
//5) The number "N", which yields the "N-th" step for the equation.
std::vector<double> eulerSolver(std::vector<std::vector<double>> &K, std::vector<std::vector<double>> &M, std::vector<double> &u_0, std::vector<double> &f, int N) {
	
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
		u_k = solveSingleStep(u_k, M, b, Pinv);
		b = createb(hf, MhK, u_k);
		std::cout << "Iternation: " << i << std::endl;
	}
	
	return u_k;
}

//Here we get the maximum step value h.
double gethmax(std::vector<std::vector<double>> &KinvM) {

	//To find h, we need to find the maximum eigenvalue of a specific matrix.
	double max_eigenvalue = findMaxEigenvalue(KinvM);
	double h_max = 2/max_eigenvalue;
	return h_max;
}

//We find the maximum eigenvalue by considering the power method.
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

//In this function we find Uk+1, having Uk. This is using the pre-conditioned gradient descent.
std::vector<double> solveSingleStep(std::vector<double> &u_k, std::vector<std::vector<double>> &A, std::vector<double> &b, std::vector<std::vector<double>> &Pinv) {

	int k = 0;
	double tol = 0.000001;
	int max_k = 1000;
	
	std::vector<double> r_k = createR(b, A, u_k);
	double dot_r0 = vecMult(r_k,r_k);
	double dot_r = 1.0;
	std::vector<double> z_k = r_k;
	double alpha = 0.0;
	while(k < max_k && dot_r > tol){
		z_k = multMatrixVec(Pinv,r_k);
		alpha = vecMult(r_k, z_k)/vecMult(z_k, multMatrixVec(A,z_k));
		generateNextStep(u_k, alpha, z_k);
		r_k = createR(b, A, u_k);
		dot_r = vecMult(r_k,r_k)/dot_r0;
		k++;
	}
	
	std::cout << "The error is " << dot_r << std::endl;
	return u_k;
}

//This function moves along the next step, based on gradient descent.
void generateNextStep(std::vector<double> &u_k, double alpha, std::vector<double> &z_k) {
	int n = u_k.size();
	for(int i = 0; i < n; i++) {
		u_k[i] += alpha*z_k[i];
	}
}

//This function find the inverse of a tri-diagonal matrix, using the LU factorisation method.
std::vector<std::vector<double>> invLUtridiagfact(std::vector<std::vector<double>> P) {
	
	int n = P.size();
	std::vector<std::vector<double>> invP = createIDmatrix(n);
	
	tdLUdecomp(P);
	tdforwardSub(P, invP);
	tdbackwardSub(P, invP);
	
	return invP;
}

//This function transforms the tri-diagonal matrix P, into strictly lower triangular and upper triangular.
void tdLUdecomp(std::vector<std::vector<double>> &P) {

	int n = P.size();
	for (int k = 0; k < n-1; k++) {
		P[k+1][k] /= P[k][k];
		P[k+1][k+1] -= P[k+1][k]*P[k][k+1]; 
	}
}

//This function performs forward substitution, so that it finds the inverse.
void tdforwardSub(std::vector<std::vector<double>> L, std::vector<std::vector<double>> &ID) {

	int n = L.size();	
	for(int columns = 0; columns < n; columns++) {
		for (int i = 1; i < n; i++) {
			ID[i][columns] -= L[i][i-1]*ID[i-1][columns];
		}
	}
}

//This function performs backward substitution, so that it finds the inverse.
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

//This function finds the inverse of a matrix using the LU factorisation.
std::vector<std::vector<double>> invLUfact(std::vector<std::vector<double>> M) {

	int n = M.size();
	std::vector<std::vector<double>> invM = createIDmatrix(n);
	
 	LUdecomp(M);
	forwardSub(M,invM);
 	backwardSub(M, invM);
	
	return invM;
}

//Finding the LU decomposition.
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

//This function performs forward substitution, so that it finds the inverse.
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

//This function performs backward substitution, so that it finds the inverse.
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
