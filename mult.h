#ifndef MULT_H
#define MULT_H

std::vector<std::vector<double>> multMatrix(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &B);

std::vector<double> multMatrixVec(std::vector<std::vector<double>> &A, std::vector<double> &v);

std::vector<double> multVecMatrix(std::vector<double> &v, std::vector<std::vector<double>> &B);

double vecMult(std::vector<double> &u, std::vector<double> &v);

#endif
