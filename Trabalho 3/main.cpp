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
void fix_bounds(std::vector<std::vector<double>>& A, bool fixed_temp = false);

// Variáveis do domínio e da simulaçaão:
constexpr double kappa {15.1};                          // coeficiente de difusividade térmica
constexpr double rho {8050};                            // massa específica
constexpr double cp {480};                              // calor específico
constexpr double g {100000};                            // geração interna

constexpr double L {0.05};                              // tamanho em x
constexpr double H {0.05};                              // tamanho em y
constexpr int Nx {101};                                 // número de nós em x
constexpr int Ny {101};                                 // número de nós em y
constexpr auto dx = L/(Nx-1);                           // comprimento da célula
constexpr auto dy = L/(Ny-1);                           // altura da célula

constexpr double T_init {20.0};                         // temperatura do meio no interior e da placa
constexpr double T_void {100.0};                        // temperatura preescrita dos contornos

constexpr double ti {0.0};                              // tempo inicial da simulação
constexpr double tf {500.0};                            // tempo final da simulação
constexpr int nsteps {4096*6};                          // número de passos de tempo
constexpr auto dt = (tf-ti)/nsteps;                     // tamanho do passo de tempo
constexpr auto gamma = g*dt/(rho*cp);
constexpr auto alfa = kappa/(rho*cp);
constexpr auto r = kappa*dt/(rho*cp*dx*dx);

int main(int argc, char* argv[]){
	std::cout << "2D NON-STEADY HEAT EQUATION SOLVER" << std::endl;

	// Variáveis auxiliares do problema:
	constexpr auto stab = (dx*dx) / (2*alfa);
	dt < stab ? std::cout << "STABILITY CRITERIA FULLFILED" << std::endl : std::cout << "STABILITY CRITERIA VIOLATED\n" << std::endl;
	std::cout << "Maximum time step allowed: " << stab << std::endl;
	std::cout << "Time step used: " << dt << std::endl; 
	std::cout << "dx = " << dx << std::endl;
	std::cout << "Coeficient r = " << r << std::endl;
	std::cout << "Total time steps to be evaluated: " << nsteps << std::endl;
	
	// Malha:
	std::vector<std::vector<double>> T(Ny, std::vector<double>(Nx, T_init));

	// Arquivo de saida:
	std::fstream save_1 {"2D_Heat.txt", std::ios::out|std::ios::trunc};

	// Inicio da iteração temporal:
	for (int step = 0; step < nsteps; step++){
		fix_bounds(T, false);
		solver_finite_difference(T, r, gamma);
		if (step%1000 == 0)
			std::cout << "Current time step: " << step << std::endl;
	}
	save_data_final(save_1, T);

	std::cout << "\nExecution reached the end";
}

void solver_finite_difference(std::vector<std::vector<double>>& T, const double r, const double gamma){
	//Calcula os valores dos nós:
	const auto N = T[0].size();
	const auto M = T.size();

	auto p = static_cast<int>(M*0.3);        // começo da área vazada: 40 % da extensão
	auto p_w = static_cast<int>(M*0.7);      // final da área vazada: 60 % da extensão
	
	// Apenas para printar uma unica vez os dados:
	static int ammount{0};
	if (ammount == 0){
		std::cout << "\nNumero de nos em x: " << N << std::endl;
		std::cout << "Numero de nos em y: " << M << std::endl;
		std::cout << "Inicio da area vazada: " << p << std::endl;
		std::cout << "Final da area vazada: " << p_w << std::endl;
	}
	ammount++;

	// Iteração nos nós do domínio
	for (int i = 1; i < N-1; i++){
		for (int j = 1; j < M-1; j++){

			// Temperatura na região interna
			if ((i>p && i<p_w) && (j>p && j<p_w)){
				continue;
			}
			else if (i==p+1 && (j>p && j<p_w)){
				T[i][j] = (1/(3 + gamma))*(gamma*T_init + 4*T[i-1][j] - T[i-2][j] );
			}
			else if (i==p_w-1 && (j>p && j<p_w)){
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

}

void save_data_final(std::fstream& u_file, const std::vector<std::vector<double>>& A){
	std::cout << "Saving data to file..." << std::endl;
	const auto N = A[0].size();
	const auto M = A.size();
	for (int i = N-1; i >= 0; i--){
		for (int j = 0; j < M; j++){
			u_file << std::setw(8) << A[i][j] << " ";
		}
		u_file << "\n";
	}
}

void fix_bounds(std::vector<std::vector<double>>& A, bool fixed_temp){

	if ((A.size() != Nx) || (A[0].size() != Ny)){
		std::cout << "Domain mismatch detected! " << std::endl;
		return;
	}

	// Arestas superior e direita:
	for (int i = 0; i < Nx; i++){
		// Superior:
		A[0][i] = T_void;
		// Direita:
		A[i][Nx-1] = T_void;
	}

	// Se a temperatura é preescrita nas arestas esquerda e inferior
	if (fixed_temp){
		for (int i = 0; i < Nx; i++){
			// Aresta esquerda:
			A[i][0] = T_void;
			// Aresta inferior:
			A[Nx-1][i] = T_void;
		}
	}
	// Se o fluxo é nulo nas arestas esquerda e inferior
	else{
		for (int i = 0; i < Nx; i++){
			// Aresta esquerda
			A[i][0] = A[i][1];
			// Aresta inferior:
			A[Nx-1][i] = A[Nx-2][i];
		}
	}

	// Garante a temperatura do fluido no interior  (entre [0.4 ~ 0.6]L x [0.4 ~ 0.6]L):
	auto p = static_cast<int>(Nx*0.3);
	auto p_w = static_cast<int>(Nx*0.7);
	for (int i = p; i < p_w; i++){
		for (int j = p; j < p_w; j++){
			A[i][j] = T_init;
		}
	}
}
