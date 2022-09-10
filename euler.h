#ifndef EULER_H
#define EULER_H

std::tuple<std::vector<double>,double,std::vector<double>> eulerSolver(std::vector<std::vector<double>> &K, std::vector<std::vector<double>> &M, std::vector<double> &u_0, std::vector<double> &f);

double gethmax(std::vector<std::vector<double>> &KinvM);

double findMaxEigenvalue(std::vector<std::vector<double>> &KinvM);

double vecMult(std::vector<double> &u, std::vector<double> &v);

std::vector<double> solveSingleStep(std::vector<double> &u_k, std::vector<std::vector<double>> &P);

std::vector<std::vector<double>> createP(std::vector<std::vector<double>> &M);

std::vector<std::vector<double>> multMatrix(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B);

std::vector<double> multMatrixVec(std::vector<std::vector<double>> &A, std::vector<double> &v);

std::vector<double> multVecMatrix(std::vector<double> &v, std::vector<std::vector<double>> &B);

std::vector<std::vector<double>> invLUfact(std::vector<std::vector<double>> M);

void LUdecomp(std::vector<std::vector<double>> &M);

void forwardSub(std::vector<std::vector<double>> L,  std::vector<std::vector<double>> &ID);

void backwardSub(std::vector<std::vector<double>> U, std::vector<std::vector<double>> &invM);

std::vector<std::vector<double>> createIDmatrix(int size);

#endif
