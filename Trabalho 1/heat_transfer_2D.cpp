#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <memory>    //std::unique_ptr
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