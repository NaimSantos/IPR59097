#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm> //std::fill
#include <fstream>   //std::fstream

#include "utilities.h"

// Protótipos de funções:
void solver_finite_difference(std::vector<std::vector<double>>& T, const double r, const double gamma);
void start_save_data(std::fstream& u_file, const std::vector<std::vector<double>>& T);
void resume_save_data(std::fstream& u_file, const std::vector<std::vector<double>>& A);
void fix_bounds(std::vector<std::vector<double>>& A, const double valueX, const double valueY);

// Variáveis do problema:
constexpr double kappa {15.1};                          // coeficiente de difusividade térmica
constexpr double rho {8050};                            // massa específica
constexpr double cp {480};                              // calor específico
constexpr double g {100000};                            // geração interna
// Variáveis do domínio:
constexpr double L {0.05};                              // tamanho em x
constexpr double H {0.05};                              // tamanho em y
constexpr int Nx {50};                                  // número de nós em x
constexpr int Ny {50};                                  // número de nós em y
constexpr auto dx = L/(Nx-1);                           // comprimento da célula
constexpr auto dy = L/(Ny-1);                           // altura da célula
// Variáveis da simulação:
constexpr double T_init {25.0};                         // temperatura inicial da placa
constexpr double ti {0.0};                              // tempo inicial da simulação
constexpr double tf {1000.0};                           // tempo final da simulação
constexpr int nsteps {4096*4};                          // número de passos de tempo
constexpr auto dt = (tf-ti)/nsteps;                     // tamanho do passo de tempo

int main(int argc, char* argv[]){
	std::cout << "2D NON-STEADY HEAT EQUATION SOLVER" << std::endl;

	// Variáveis auxiliares do problema:
	constexpr auto gamma = g*dt/(rho*cp);
	constexpr auto alfa = kappa/(rho*cp);
	constexpr auto r = kappa*dt/(rho*cp*dx*dx);
	constexpr auto stab = (dx*dx) / (4*alfa);
	dt < stab ? std::cout << "Estabilidade satisfeista" << std::endl : std::cout << "Estabilidade violada" << std::endl;
	std::cout << "dt\t\tstability" << std::endl; 
	std::cout << dt << "\t" << stab << std::endl; 
	std::cout << "\ndx = " << dx << "\nr = " << r << "\ngamma = " << gamma << "\nalfa = " << alfa << std::endl;

	// Malha:
	std::vector<std::vector<double>> T(Ny, std::vector<double>(Nx, T_init));
	fix_bounds(T, 100, 100);
	// Arquivo de saida:
	std::fstream save_1 {"2D_Heat.dat", std::ios::out|std::ios::trunc};
	std::fstream save_2 {"2D_Heat.dat", std::ios::out|std::ios::app};

	// Inicio da execução:
	start_save_data(save_1, T);
	save_1.close(); //Fclose
	for (int step = 0; step < nsteps; step++){
		solver_finite_difference(T, r, gamma);
		fix_bounds(T, 100, 100);
		if (step%4096==0)
			resume_save_data(save_2, T);
	}

	std::cout << "\nExecution reached the end" << std::endl;
}

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
	u_file << "\n";

	const auto N = A[0].size();
	const auto M = A.size();
	for (int i = N-1; i >= 0; i--){
		for (int j = 0; j < M; j++){
			u_file << std::setw(8) << A[i][j] << " ";
		}
		u_file << "\n";
	}
}

void fix_bounds(std::vector<std::vector<double>>& A, const double valueX, const double valueY){
	auto m = A.size();
	auto n = A[0].size();

	// Aresta superior:
	std::fill(A[0].begin(), A[0].end(), valueX);

	for (int k = 0; k < m; k++){
		// Aresta direita
		A[k][0] = A[k][1];
		// Aresta esquerda:
		A[k][n-1] = valueY;
		// Aresta inferior:
		A[m-1][k] = A[m-2][k];
	}
}
