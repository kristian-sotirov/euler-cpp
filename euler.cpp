#include <iostream>
#include <vector>
#include <tuple>
#include "euler.h"
using namespace std;

//THIS FUNCTION IS UNFINISHED!!!!!!!!!!
tuple<vector<double>, double, vector<double>> eulerSolver(vector<vector<double>> &K, vector<vector<double>> &M, vector<double> &u_0, vector<double> &f) {
	
	auto return_value = make_tuple(u_0, 0.0, f);
	return return_value;

}

double gethmax(vector<vector<double>> KinvM) {
	max_eigenvalue = findMaxEigenvalue(KinvM);
	h_max = 2/max_eigenvalue;
	return h_max;
}

//THIS FUNCTION IS UNFINISHED!!!!!!!!!
double findMaxEigenvalue(KinvM) {
	return 1.0;
}
//THIS FUNCTION IS UNFINISHED!!!!!!!!!
vector<double> solveSingleStep(vector<double> &u_k vector<vector<double>> &P) {
	vector<double> solution_vector(5, 0.0);
	return solution_vector;
}
//THIS FUNCTION IS UNFINISHED!!!!!!!!!
vector<vector<double>> createP(vector<vector<double>> &M) {
	vector<vector<double>> P = {{0.0,0.0,0.0}, {0.0,0.0,0.0}}
	return P;
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
