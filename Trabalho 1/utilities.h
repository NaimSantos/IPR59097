#pragma once

#include <vector>
#include <cmath>
#include <iostream>

template<typename T>
bool is_diagonal_dom(const std::vector<std::vector<T>>& A){
	auto n_rows = A.size();
	auto n_col = A[0].size();
	if (n_rows != n_col)
		return false;

	for (int i = 0; i < n_rows; i++){
		double sum = 0.0;
		double diag = std::fabs(A[i][i]);
		for (int j = 0; j < n_col; j++){
			if (i == j)
				continue;
			sum += std::fabs(A[i][j]);
		}
		if (sum > diag)
			return false;
	}

	return true;
}

template<typename T>
void GS_Solver(const std::vector<std::vector<T>>& A, std::vector<T>& B, std::vector<T>& X){
	// Caso não seja diagonal dominante, a convergência não é garantida
	if (!is_diagonal_dom(A))
		return;

	// Dimensões não são compatíveis:
	if ((A.size() != A[0].size()) || A.size() != B.size())
		return;

	auto Y = B;
	auto E = X;

	int counter {0};
	bool teste = false;
	const double eps {0.000001};

	int m = A.size();
	int n = A[0].size();
	while(!teste && counter<40){
		std::cout << "Iter = " << counter << std::endl;
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

template <typename T>
void print_array_2D(const std::vector<std::vector<T>> M){
	auto nrow = M.size();
	auto ncol = M[0].size();
	for(int i = 0; i < nrow; i++){
		for(int j = 0; j < ncol; j++){
			std::cout << M[i][j] << ' ';
		}
		std::cout << '\n';
	}
}

template <typename T>
void print_array_1D(const std::vector<T> M){
	auto nrow = M.size();
	for(int i = 0; i < nrow; i++){
		std::cout << M[i] << ' ';
	}
	std::cout << '\n';
}
