#ifndef EULER_H
#define EULER_H

tuple<std::vector<int>,int,std::vector<int>> eulerSolver(vector<vector<double>> &K, vector<vector<double>> &M, vector<double> &u_0, vector<double> &f);

double gethmax(vector<vector<double>> KinvM);

double findMaxEigenvalue(KinvM);

vector<double> solveSingleStep(vector<double> &u_k vector<vector<double>> &P);

vector<vector<double>> createP(vector<vector<double>> &M);

vector<vector<double>> multMatrix(vector<vector<double>> &A, vector<vector<double>> &B);


#endif
