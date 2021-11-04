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
void save_data_final(std::fstream& u_file, const std::vector<std::vector<double>>& T);
void fix_bounds(std::vector<std::vector<double>>& T, bool dirichlet = false);
void print_messages();
void try_static_assertions();

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
constexpr double T_init {20.0};                         // temperatura inicial da placa
constexpr double T_fixed {100.0};                       // temperatura prescrita no contorno
constexpr double T_inf {20.0};                          // temperatura do fluido no interior
constexpr double ti {0.0};                              // tempo inicial da simulação
constexpr double tf {500.0};                            // tempo final da simulação
constexpr auto alfa = kappa/(rho*cp);
constexpr auto stab = (dx*dx) / (2*alfa);               // máximo passo de tempo para estabilidade
constexpr auto dt = 0.95*stab;                          // passo de tempo (95 % do passo máximo)
constexpr int nsteps = 1 + static_cast<int>((tf-ti)/dt);// número de passos de tempo
constexpr auto gamma = g*dt/(rho*cp);
constexpr auto r = kappa*dt/(rho*cp*dx*dx);

int main(int argc, char* argv[]){

	try_static_assertions();
	print_messages();

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

	std::cout << "\n------- Execution reached the end.------- ";
}


void solver_finite_difference(std::vector<std::vector<double>>& T, const double r, const double gamma){

	auto p = static_cast<int>(Nx*0.3);        // começo da área vazada: 30 % da extensão
	auto p_w = static_cast<int>(Nx*0.7);      // final da área vazada: 70 % da extensão

	// Apenas para printar uma unica vez os dados:
	static int ammount{0};
	if (ammount == 0){
		std::cout << "\nNumero de nos em x: " << Nx << std::endl;
		std::cout << "Numero de nos em y: " << Ny << std::endl;
		std::cout << "Inicio da area vazada: " << p << std::endl;
		std::cout << "Final da area vazada: " << p_w << std::endl;
	}
	ammount++;

	// Iteração nos nós do domínio, exceto contornos externos
	for (int i = 1; i < Nx-1; i++){
		for (int j = 1; j < Ny-1; j++){
			// Temperatura na região vazada
			if ((i>p && i<p_w) && (j>p && j<p_w)){
				continue;
			}
			// Temperatura nas regiões do contorno interno (área vazada)
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
			// Temperatura nos demais nós
			else{
				T[i][j] = T[i][j] + r*(T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1] - 4*T[i][j]) + gamma;
			}
		}
	}

}

void fix_bounds(std::vector<std::vector<double>>& T, bool dirichlet){

	// Arestas superior e direita:
	for (int i = 0; i < Nx; i++){
		// Superior:
		T[0][i] = T_fixed;
		// Direita:
		T[i][Nx-1] = T_fixed;
	}

	// Temperatura prescrita nos contornos esquerdo e inferior (Dirichlet)
	if (dirichlet){
		for (int i = 0; i < Nx; i++){
			// Aresta esquerda:
			T[i][0] = T_fixed;
			// Aresta inferior:
			T[Nx-1][i] = T_fixed;
		}
	}
	// Fluxo nulo nos contornos esquerdo e inferior (Neumann)
	else{
		for (int i = 0; i < Nx; i++){
			// Aresta esquerda
			T[i][0] = T[i][1];
			// Aresta inferior:
			T[Nx-1][i] = T[Nx-2][i];
		}
	}

	// Garante a temperatura do fluido no interior  (entre [0.3 ~ 0.7]L x [0.3 ~ 0.7]L):
	auto p = static_cast<int>(Nx*0.3);
	auto p_w = static_cast<int>(Nx*0.7);
	for (int i = p; i < p_w; i++){
		for (int j = p; j < p_w; j++){
			T[i][j] = T_inf;
		}
	}
}

void save_data_final(std::fstream& u_file, const std::vector<std::vector<double>>& T){

	std::cout << "Saving data to file..." << std::endl;
	const auto N = T[0].size();
	const auto M = T.size();
	for (int i = N-1; i >= 0; i--){
		for (int j = 0; j < M; j++){
			u_file << std::setw(8) << T[i][j] << " ";
		}
		u_file << "\n";
	}
}

void print_messages(){

	std::cout << "------- 2D NON-STEADY HEAT EQUATION SOLVER -------" << std::endl;
	dt < stab ? std::cout << "\nSTABILITY CRITERIA FULLFILED!" << std::endl : std::cout << "\nSTABILITY CRITERIA VIOLATED!" << std::endl;
	std::cout << "Maximum time step allowed: " << stab << std::endl;
	std::cout << "Time step used: " << dt << std::endl; 
	std::cout << "\ndx = " << dx << std::endl;
	std::cout << "Coeficient r = " << r << std::endl;
	std::cout << "Total time steps to be evaluated: " << nsteps << std::endl;
}

void try_static_assertions(){

	static_assert(Nx != 0 && Ny != 0, "Numero de nos nao pode ser nulo");
	static_assert(Nx == Ny, "Numero de nos precisa em x e y ser igual");
	static_assert(Nx > 10 || Ny > 10, "Numero de nos precisa ser maior que 10");
}
