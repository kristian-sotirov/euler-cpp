#ifndef EULER_H
#define EULER_H

std::tuple<std::vector<double>,double,std::vector<double>> eulerSolver(std::vector<std::vector<double>> &K, std::vector<std::vector<double>> &M, std::vector<double> &u_0, std::vector<double> &f);

double gethmax(std::vector<std::vector<double>> &KinvM);

double findMaxEigenvalue(std::vector<std::vector<double>> &KinvM);

std::vector<double> solveSingleStep(std::vector<double> &u_k, std::vector<std::vector<double>> &P);

std::vector<std::vector<double>> createP(std::vector<std::vector<double>> &M);

std::vector<std::vector<double>> multMatrix(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B);


#endif
