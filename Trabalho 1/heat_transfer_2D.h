#pragma once

void start_save_data(std::fstream& u_file, const std::vector<std::vector<double>>& T);
void solver_finite_difference(std::vector<std::vector<double>>& T, const double r, const double gamma);
void set_bounds(std::vector<std::vector<double>>& A, double value);
void set_full_bound(std::vector<std::vector<double>>& A, double value, int line);
void resume_save_data(std::fstream& u_file, const std::vector<std::vector<double>>&);
