#include <iostream>
#include <vector>
#include "euler.h"
using namespace std;


tuple<vector<double>, double, vector<double>> eulerSolver(vector<vector<double>> &K, vector<vector<double>> &M, vector<double> &u_0, vector<double> &f) {



}

double gethmax(vector<vector<double>> KinvM) {

}

double findMaxEigenvalue(KinvM) {

}

vector<double> solveSingleStep(vector<double> &u_k vector<vector<double>> &P) {

}

vector<vector<double>> createP(vector<vector<double>> &M) {

}

vector<vector<double>> multMatrix(vector<vector<double>> &A, vector<vector<double>> &B) {

	int m = A.size();
	int p = A[0].size();
	
	if(p != B.size()) {
		cout << "The matrices cannot be multiplied";
		return A;
	}
	int n = B[0].size();
	
	vector<vector<double>> C(m, vector<double>(n));
	
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < n; j++) {
		
			for(int k = 0; k < p; k++) {
				
				C[i][j] += А[i][k]*B[k][j];
			}
		}
	}
	
	return C;
}
