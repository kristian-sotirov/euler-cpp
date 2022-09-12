#ifndef MATRIXGENERATION_H
#define MATRIXGENERATION_H

std::vector<double> createR(std::vector<double> &b, std::vector<std::vector<double>> &A, std::vector<double> &x);

std::vector<std::vector<double>> createP(std::vector<std::vector<double>> &M);

std::vector<std::vector<double>> createMhK(std::vector<std::vector<double>> &M, double h, std::vector<std::vector<double>> &K);

std::vector<double> createhf(double h, std::vector<double> &f);

std::vector<double> createb(std::vector<double> &hf, std::vector<std::vector<double>> &MhK, std::vector<double> &u_k);

#endif
