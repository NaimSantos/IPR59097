#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>   //std::fstream

#include "utilities.h"
#include "heat_transfer_2D.h"

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

	// Variáveis auxiliares do problema:
	constexpr auto gamma = g*dt/(rho*cp);
	constexpr auto alfa = kappa/(rho*cp);
	constexpr auto r = kappa*dt/(rho*cp*dx*dx);
	constexpr auto stab = (dx*dx) / (4*alfa);
	std::cout << "stability = " << stab << "\ndt = "<< dt << ", dx = " << dx << ", r = " << r << ", gamma = " << gamma << ", alfa = " << alfa << std::endl;


	// Malha:
	std::vector<std::vector<double>> T(Ny, std::vector<double>(Nx, T_init));
	fix_bounds(T, 100, 100);
	// Arquivo de saida:
	std::fstream save_1 {"2D_Heat.dat", std::ios::out|std::ios::trunc};
	std::fstream save_2 {"2D_Heat.dat", std::ios::out|std::ios::app};

	// Inicio da execução:
	start_save_data(save_1, T);

	for (int step = 0; step < nsteps; step++){
		solver_finite_difference(T, r, gamma);
		fix_bounds(T, 100, 100);
		if (step%4096==0)
			resume_save_data(save_2, T);
	}


	std::cout << "\nExecution reached the end" << std::endl;
}
