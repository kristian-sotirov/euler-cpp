#include <iostream>
#include <vector>
#include <cmath>
#include "mult.h"

double vecMult(std::vector<double> u, std::vector<double> v) {
	
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

std::vector<double> multMatrixVec(std::vector<std::vector<double>> &A, std::vector<double> v) {
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
