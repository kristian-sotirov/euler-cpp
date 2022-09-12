#ifndef EULER_H
#define EULER_H

std::tuple<std::vector<double>,double,std::vector<int>> eulerSolver(std::vector<std::vector<double>> &K, std::vector<std::vector<double>> &M, std::vector<double> &u_0, std::vector<double> &f, int N);

double gethmax(std::vector<std::vector<double>> &KinvM);

double findMaxEigenvalue(std::vector<std::vector<double>> &KinvM);

std::tuple<std::vector<double>, int> solveSingleStep(std::vector<double> &u_k, std::vector<std::vector<double>> &A, std::vector<double> &b, std::vector<std::vector<double>> &P);

std::vector<std::vector<double>> createP(std::vector<std::vector<double>> &M);

std::vector<std::vector<double>> createMhK(std::vector<std::vector<double>> &M, double h, std::vector<std::vector<double>> &K);

std::vector<double> createhf(double h, std::vector<double> &f);

std::vector<double> createb(std::vector<double> &hf, std::vector<std::vector<double>> &MhK, std::vector<double> &u_k);

std::vector<std::vector<double>> invLUtridiagfact(std::vector<std::vector<double>> P);

void tdLUdecomp(std::vector<std::vector<double>> &P);

void tdforwardSub(std::vector<std::vector<double>> L, std::vector<std::vector<double>> &ID);

void tdbackwardSub(std::vector<std::vector<double>> U, std::vector<std::vector<double>> &invP);

std::vector<std::vector<double>> invLUfact(std::vector<std::vector<double>> M);

void LUdecomp(std::vector<std::vector<double>> &M);

void forwardSub(std::vector<std::vector<double>> L,  std::vector<std::vector<double>> &ID);

void backwardSub(std::vector<std::vector<double>> U, std::vector<std::vector<double>> &invM);

std::vector<std::vector<double>> createIDmatrix(int size);

#endif
