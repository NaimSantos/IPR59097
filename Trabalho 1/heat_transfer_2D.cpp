#include <iostream>
#include <cmath>
#include <string>
#include <algorithm> //std::fill
#include <fstream>   //std::fstream

#include "utilities.h"

void solver_finite_difference(std::vector<std::vector<double>>& T, const double r, const double gamma){
	//Calcula os valores dos nós:
	const auto N = A[0].size();
	const auto M = A.size();
	for (int i = 1; i < N-1; i++){
		for(int j = 1; j < M -1; j++){
			T[i][j] = r*(T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1] - 4*T[i][j]) + T[i1][j] + gamma;
		}
	}
}

void start_save_data(std::fstream& u_file, const std::vector<std::vector<double>>& T){
	u_file << "Perfil de Temperatura\n";
	u_file << "x y T\n";
	const auto N = A[0].size();
	const auto M = A.size();
	for (int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			//if (j == M-1)
			//	u_file << T[i][j];
			//else
				u_file << T[i][j] << " ";
		}
		u_file << "\n\n";
	}

}

void resume_save_data(std::fstream& u_file, const std::vector<std::vector<double>>&){
	const auto N = A[0].size();
	const auto M = A.size();
	for (int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
				u_file << T[i][j] << " ";
		}
		u_file << "\n\n";
	}
}

void fix_bounds_newman(std::vector<std::vector<double>>& A){
	auto m = A.size();
	auto n = A[0].size();
	
}

void fix_bounds_dirichillet(std::vector<std::vector<double>>& A, const double valueX, const double valueY){
	auto m = A.size();
	auto n = A[0].size();
	std::fill(A[0].begin(), A[0].end(), valueX);
	for (int i = 0; i < m; i++)
		A[m-1][i] = valueY;
}
