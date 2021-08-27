#pragma once

void start_save_data(std::fstream& printer);
void solver_finite_difference(std::vector<std::vector<double>>&);
void set_bounds(std::vector<std::vector<double>>& A, double value);
void set_full_bound(std::vector<std::vector<double>>& A, double value, int line);
void save_data(const std::vector<std::vector<double>>&, std::fstream& filename);
