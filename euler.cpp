#include <iostream>
#include <vector>
#include <tuple>
#include "euler.h"
#include <string>

//THIS FUNCTION IS UNFINISHED!!!!!!!!!!
std::tuple<std::vector<double>, double, std::vector<double>> eulerSolver(std::vector<std::vector<double>> &K, std::vector<std::vector<double>> &M, std::vector<double> &u_0, std::vector<double> &f) {
	
	std::tuple<std::vector<double>, double, std::vector<double>> return_value = std::make_tuple(u_0, 0.0, f);
	
	std::vector<std::vector<double>> invM = invLUfact(M);
	
	return return_value;

}

double gethmax(std::vector<std::vector<double>> &KinvM) {
	double max_eigenvalue = findMaxEigenvalue(KinvM);
	double h_max = 2/max_eigenvalue;
	return h_max;
}

//THIS FUNCTION IS UNFINISHED!!!!!!!!!
double findMaxEigenvalue(std::vector<std::vector<double>> &KinvM) {
	return 1.0;
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


std::vector<std::vector<double>> invLUfact(std::vector<std::vector<double>> &M) {
	
	// Done
	std::vector<std::vector<double>> M_decompLU = LUdecomp(M);
	
	// Done
	std::vector<std::vector<double>> semi_invM = forwardSub(M_decompLU);
	
	
	std::vector<std::vector<double>> invM = backwardSub(M_decompLU, semi_invM);
	
	return invM;
}


std::vector<std::vector<double>> LUdecomp(std::vector<std::vector<double>> &M) {
	std::vector<std::vector<double>> M_decompLU = M;
	int n = M.size();
	for (int k = 0; k < n-1; k++) {
		for (int j = k+1; j < n; j++) {
			M_decompLU[j][k] /= M_decompLU[k][k];
			for (int i = k+1; i < n; i++) {
				M_decompLU[j][i] -= M_decompLU[j][k]*M_decompLU[k][i];
			}
		}
	}
	return M_decompLU;
}

std::vector<std::vector<double>> forwardSub(std::vector<std::vector<double>> &L) {
	int n = L.size();
	std::vector<std::vector<double>> ID = createIDmatrix(n);
	//std::vector<std::vector<double>> semi_invM = std::vector<std::vector<double>>(size, std::vector<double>(size, 0.0));
	
	for(int IDcolumn = 0; IDcolumn < n; IDcolumn++) {
		//semi_invM[0][IDcolumn] = ID[0][IDcolumn];
		for (int i = 1; i < n; i++) {
			for (int j = 0; j < i; j++) {
				ID[i][IDcolumn] -= L[i][j]*ID[j][IDcolumn];
			}
			//ID[i][IDcolumn] /= L[i][i];
		}
	}
	
	return ID;
}

//THIS FUNCTION IS UNFINISHED!!!!!!!!!
std::vector<std::vector<double>> backwardSub(std::vector<std::vector<double>> &U, std::vector<std::vector<double>> &invM) {
	
	return invM;
}


std::vector<std::vector<double>> createIDmatrix(int size) {

	std::vector<std::vector<double>> ID = std::vector<std::vector<double>>(size, std::vector<double>(size, 0.0));
	for(int i = 0; i < size; i++) {
		ID[i][i] = 1.0;
	}
	return ID;
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
