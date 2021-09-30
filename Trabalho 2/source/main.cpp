#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <chrono>

void GS_Solver(const std::vector<std::vector<double>>& A, std::vector<double>& B, std::vector<double>& X);
template <typename T> void linspace(std::vector<double>& V, int Num, T xf, T xi);
void save_data_2D(const std::vector<std::vector<double>>& A, const std::string& filename);
template <typename T>
void print_array_2D(const std::vector<std::vector<T>> M);

constexpr double kappa {0.6};                            // coeficiente de condutividade térmica
constexpr double h {100.0};                               // coeficiente de troca de calor por convecção
constexpr double cp {1200.0};                            // calor específico
constexpr double rho {600.0};                            // massa específica
constexpr double g {100000.0};                           // geração interna
constexpr double T_inf{20.0};                            // temperatura do ambiente
constexpr double T_0 {20};                               // temperatura inicial da placa
constexpr double L {0.03};                               // comprimento da placa
constexpr int N {11};                                    // número de nós em x
constexpr double ti {0.0};                               // tempo inicial
constexpr double tf {5.0};                              // tempo total de simulação
constexpr double dt {1};                               // passo de tempo
constexpr auto dx = L /(N - 1);                          // intervalo
constexpr auto nsteps = static_cast<int>((tf - ti)/dt);  // número de passos de tempo

int main(int argc, char* argv[]){
	constexpr auto alpha = kappa / (rho * cp);
	constexpr auto r = (alpha * dt)/(dx*dx);
	constexpr auto gamma = (alpha * g * dt) / kappa;
	constexpr auto mu = (2 * dx * h) / kappa;

	std::cout << "alpha = " << alpha << std::endl;
	std::cout << "r = " << r << std::endl;
	std::cout << "gamma = " << gamma << std::endl;
	std::cout << "mi = " << mu << std::endl;

	std::vector<std::vector<double>> T(nsteps, std::vector<double>(N, T_0));  // matriz N x nsteps, para as temperaturas em todos os tempos
	std::vector<double> Pos(N, 0.0);                                          // vetor com as posições avaliadas
	linspace(Pos, N, L, 0.0);

	print_array_2D(T);
	for (int i = 1; i < nsteps; i++){
		T[i][0] = -r*T[i][0] + 2*r*T[i][1] + gamma;
		for (int j = 1; j < N-1; j++){
			
			T[i][j] = r*T[i][j-1] - r*T[i][j] + r*T[i][j+1] + gamma;
			std::cout << "i = " << i << " j = " << j << "\tT = " << T[i][j] << std::endl;
		}
		T[i][N-1] = 2*r*T[i][N-2] - (r + mu)*T[i][N-1] + r*mu*T_inf + gamma;
	}
	const std::string file1 {"output_data.txt"};
	save_data_2D(T, file1);
	print_array_2D(T);
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
void linspace(std::vector<double>& V, int Num, T xf, T xi){
	auto h = (xf - xi) / (Num-1);
	auto n = static_cast<int>(V.size());
	for (int i = 0; i < n; i++){
		V[i] = xi + i*h;
	}
}

void save_data_2D(const std::vector<std::vector<double>>& A, const std::string& filename){
	std::fstream saver {filename, std::ios::out|std::ios::trunc};
	int n = A[0].size();
	int m = A.size();
	
	for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
			saver << A[i][j] << ' ';
		}
		saver << '\n';
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