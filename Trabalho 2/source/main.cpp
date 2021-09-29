#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <chrono>

void GS_Solver(const std::vector<std::vector<T>>& A, std::vector<T>& B, std::vector<T>& X);
template <typename T>
void linspace(std::vector<double>& V, const int Num, const T xf, const T xi);

int main(int argc, char* argv[]){


	std::cout << "Execution reached the end" << std::endl;
}

void GS_Solver(const std::vector<std::vector<double>>& A, std::vector<double>& B, std::vector<double>& X){
	// Dimensões não são compatíveis:
	if ((A.size() != A[0].size()) || A.size() != B.size())
		return;

	auto Y = B;
	auto E = X;

	int counter {0};
	bool teste = false;
	const double eps {0.0000001};

	int m = A.size();
	int n = A[0].size();
	while(!teste && counter<(m*5)){
		teste = true;
		for (int i = 0; i < m; i++){
			Y[i] = (B[i] / A[i][i]);
			for (int j = 0; j < n; j++){
				if (j == i)
					continue;
				Y[i] = Y[i] - ((A[i][j] / A[i][i]) * X[j]);
				X[i] = Y[i]; // Escreve em X a estimativa encontrada
			}
			auto res = std::fabs(((X[i] - E[i]) / X[i])) <= eps;
			teste = teste & res;
			E[i] = X[i];
		}
		counter++;
	}
}

template <typename T>
void linspace(std::vector<double>& V, const int Num, const T xf, const T xi){
	auto h = (xf - xi) / (Num-1);
	auto n = static_cast<int>(V.size());
	for (int i = 0; i < n; i++){
		V[i] = xi + i*h;
	}
}
