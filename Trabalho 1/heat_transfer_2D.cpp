#include <iostream>
#include <cmath>
#include <string>
#include <algorithm> //std::fill
#include <fstream>   //std::fstream

#include "utilities.h"

void solver_finite_difference(std::vector<std::vector<double>>&){
	auto pi_teste = 3.14;
	pi_teste = std::sqrt(pi_teste);
}

void start_save_data(std::fstream& printer){
	std::cout << "ECho" << std::endl;
	printer << "Perfil de Temperatura\n";
	printer << "x y T\n";
}

void save_data(const std::vector<std::vector<double>>&, std::fstream& filename){
	std::cout << "Execução atingindo a funcao save data" << std::endl;
}

void fix_bounds(std::vector<std::vector<double>>& A, const double value){
	auto m = A.size();
	auto n = A[0].size();
	A[0][0] = value;
	A[m][0] = value;
	A[0][n] = value;
	A[m][n] = value;
}
