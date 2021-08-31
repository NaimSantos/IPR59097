#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <algorithm> //std::fill
#include <fstream>   //std::fstream

#include "utilities.h"

void solver_finite_difference(std::vector<std::vector<double>>& T, const double r, const double gamma){
	//Calcula os valores dos nós:
	const auto N = T[0].size();
	const auto M = T.size();
	for (int i = 1; i < N-1; i++){
		for (int j = 1; j < M-1; j++){
			T[i][j] = r*(T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1] - 4*T[i][j]) + T[i][j] + gamma;
		}
	}
}

void start_save_data(std::fstream& u_file, const std::vector<std::vector<double>>& A){
	u_file << "Perfil de Temperatura\n";
	
	const auto N = A[0].size();
	const auto M = A.size();
	for (int i = 0; i < N; i++){
		for (int j = 0; j < M; j++){
			u_file << std::setprecision(3) << std::setw(8) << A[i][j] << " ";
		}
		u_file << "\n";
	}
}

void resume_save_data(std::fstream& u_file, const std::vector<std::vector<double>>& A){
	static int count = 0;
	const auto N = A[0].size();
	const auto M = A.size();
	u_file << "\n";
	for (int i = N-1; i >= 0; i--){
		for (int j = 0; j < M; j++){
			u_file << std::setw(8) << A[i][j] << " ";
		}
		u_file << "\n";
	}
	count++;
}

void fix_bounds(std::vector<std::vector<double>>& A, const double valueX, const double valueY){
	auto m = A.size();
	auto n = A[0].size();

	// Aresta superior:
	std::fill(A[0].begin(), A[0].end(), valueX); 
	// Aresta direita
	for (int k = 0; k < m; k++)
		A[k][0] = A[k][1];
	// Aresta esquerda:
	for (int k = 0; k < m; k++)
		A[k][n-1] = valueY;
	// Aresta inferior:
	for (int k = 0; k < m; k++)
		A[m-1][k] = A[m-2][k];
}
