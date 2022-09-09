#include <iostream>
#include <vector>
#include <tuple>
#include "euler.h"
#include <string>

//THIS FUNCTION IS UNFINISHED!!!!!!!!!!
std::tuple<std::vector<double>, double, std::vector<double>> eulerSolver(std::vector<std::vector<double>> &K, std::vector<std::vector<double>> &M, std::vector<double> &u_0, std::vector<double> &f) {
	
	std::tuple<std::vector<double>, double, std::vector<double>> return_value = std::make_tuple(u_0, 0.0, f);
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
