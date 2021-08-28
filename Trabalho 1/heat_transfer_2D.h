#pragma once

void solver_finite_difference(std::vector<std::vector<double>>& T, const double r, const double gamma);
void start_save_data(std::fstream& u_file, const std::vector<std::vector<double>>& T);
void resume_save_data(std::fstream& u_file, const std::vector<std::vector<double>>& A);
void fix_bounds(std::vector<std::vector<double>>& A, const double valueX, const double valueY);