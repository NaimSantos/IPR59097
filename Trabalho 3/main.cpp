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
void save_data_final(std::fstream& u_file, const std::vector<std::vector<double>>& A);
void fix_bounds(std::vector<std::vector<double>>& A, const double valueT, bool fixed_temp = false);

// Variáveis do problema:
constexpr double kappa {15.1};                          // coeficiente de difusividade térmica
constexpr double rho {8050};                            // massa específica
constexpr double cp {480};                              // calor específico
constexpr double g {100000};                            // geração interna
// Variáveis do domínio:
constexpr double L {0.05};                              // tamanho em x
constexpr double H {0.05};                              // tamanho em y
constexpr int Nx {100};                                  // número de nós em x
constexpr int Ny {100};                                  // número de nós em y
constexpr auto dx = L/(Nx-1);                           // comprimento da célula
constexpr auto dy = L/(Ny-1);                           // altura da célula
// Variáveis da simulação:
constexpr double T_init {20.0};                         // temperatura inicial da placa
constexpr double T_void {100.0};                         // temperatura inicial da placa

constexpr double ti {0.0};                              // tempo inicial da simulação
constexpr double tf {1000.0};                           // tempo final da simulação
constexpr int nsteps {4096*8};                          // número de passos de tempo
constexpr auto dt = (tf-ti)/nsteps;                     // tamanho do passo de tempo
constexpr auto gamma = g*dt/(rho*cp);
constexpr auto alfa = kappa/(rho*cp);
constexpr auto r = kappa*dt/(rho*cp*dx*dx);

int main(int argc, char* argv[]){
	std::cout << "2D NON-STEADY HEAT EQUATION SOLVER" << std::endl;

	// Variáveis auxiliares do problema:

	constexpr auto stab = (dx*dx) / (2*alfa);
	dt < stab ? std::cout << "Estabilidade satisfeista" << std::endl : std::cout << "Estabilidade violada" << std::endl;
	std::cout << "dt\t\tstability" << std::endl; 
	std::cout << dt << "\t" << stab << std::endl; 
	std::cout << "\ndx = " << dx << "\nr = " << r << std::endl;

	// Malha:
	std::vector<std::vector<double>> T(Ny, std::vector<double>(Nx, T_init));
	fix_bounds(T, T_void, true);
	// Arquivo de saida:
	std::fstream save_1 {"2D_Heat.txt", std::ios::out|std::ios::trunc};

	// Inicio da execução:

	for (int step = 0; step < nsteps; step++){
		solver_finite_difference(T, r, gamma);
		fix_bounds(T, T_void, true);	
	}
	save_data_final(save_1, T);

	std::cout << "\nExecution reached the end";
}

void solver_finite_difference(std::vector<std::vector<double>>& T, const double r, const double gamma){
	//Calcula os valores dos nós:
	const auto N = T[0].size();
	const auto M = T.size();
	static int ammount{0};
	auto p = static_cast<int>(M*0.4);
	auto p_w = static_cast<int>(M*0.6);
	if (ammount == 0){
		std::cout << "N: " << N << std::endl;
		std::cout << "M: " << M << std::endl;
		std::cout << "p: " << p << std::endl;
		std::cout << "p_w: " << p_w << std::endl;
	}
	for (int i = 1; i < N-1; i++){
		for (int j = 1; j < M-1; j++){

			// Temperatura na região interna
			if ((i>p && i<p_w) && (j>p && j<p_w)){
				continue;
			}
			else if (i==p+1 && (j>=p && j<=p_w)){
				T[i][j] = (1/(3 + gamma))*(gamma*T_init + 4*T[i-1][j] - T[i-2][j] );
			}
			else if (i==p_w-1 && (j>=p && j<=p_w)){
				T[i][j] = (1/(gamma - 3))*(gamma*T_init - 4*T[i+1][j] + T[i+2][j] );
			}
			else if (j==p+1 && (i>=p && i<=p_w)){
				T[i][j] = (1/(3 + gamma))*(gamma * T_init + 4*T[i][j-1] - T[i][j-2] );
			}
			else if(j==p_w-1 && (i>=p && i<=p_w)){
				T[i][j] = (1/(gamma - 3))*(gamma*T_init - 4*T[i][j+1] + T[i][j+2]) ;
			}
			else{
				T[i][j] = T[i][j] + r*(T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1] - 4*T[i][j]) + gamma;
			}
		}
	}
	ammount++;
}


void save_data_final(std::fstream& u_file, const std::vector<std::vector<double>>& A){
	const auto N = A[0].size();
	const auto M = A.size();
	for (int i = N-1; i >= 0; i--){
		for (int j = 0; j < M; j++){
			u_file << std::setw(8) << A[i][j] << " ";
		}
		u_file << "\n";
	}
}

void fix_bounds(std::vector<std::vector<double>>& A, const double valueT, bool fixed_temp){
	auto m = A.size();
	auto n = A[0].size();
	
	auto p = static_cast<int>(m*0.4);
	auto p_w = static_cast<int>(m*0.6);


	// Aresta superior:
	std::fill(A[0].begin(), A[0].end(), valueT);
	
	if (fixed_temp){
		for (int k = 0; k < m; k++){
			// Aresta direita
			A[k][0] = T_void;
			// Aresta esquerda:
			A[k][n-1] = valueT;
			// Aresta inferior:
			A[m-1][k] = T_void;
		}
	}
	else{
		for (int k = 0; k < m; k++){
			// Aresta direita
			A[k][0] = A[k][1];
			// Aresta esquerda:
			A[k][n-1] = valueT;
			// Aresta inferior:
			A[m-1][k] = A[m-2][k];
		}
	}



	for (int i = p; i < p_w; i++){
		for (int j = p; j < p_w; j++){
			A[i][j] = T_init;
		}
	}
	
}
