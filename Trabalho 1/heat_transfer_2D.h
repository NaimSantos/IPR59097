#pragma once

//Protótipos das funções:
void start_save_data(std::fstream& printer);
void solver_finite_difference(std::vector<std::vector<double>>&);
//void fix_bounds(std::vector<std::vector<double>>&, const double value);
void save_data(const std::vector<std::vector<double>>&, std::fstream& filename);