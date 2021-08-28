#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>   //std::fstream

#include "utilities.h"
#include "heat_transfer_2D.h"

// Variáveis do problema:
constexpr double kappa {0.6};                          // coeficiente de difusividade térmica
constexpr double rho {1000};                           // massa específica
constexpr double cp {0.6};                             // calor específico
// Variáveis do domínio:
constexpr double L {1.0};                              // tamanho em x
constexpr double H {1.0};                              // tamanho em y
constexpr int Nx {10};                                 // número de nós em x
constexpr int Ny {10};                                 // número de nós em y
constexpr auto dx = L/(Nx-1);                          // comprimento da célula
constexpr auto dy = L/(Ny-1);                          // altura da célula
// Variáveis da simulação:
constexpr double T_init {20.0};                        // temperatura inicial da placa
constexpr double ti {0.0};                             // tempo inicial da simulação
constexpr double tf {100.0};                           // tempo final da simulação
constexpr int nsteps {4096};                           // número de passos de tempo
constexpr auto dt = (tf-ti)/nsteps;                    // tamanho do passo de tempo

int main(int argc, char* argv[]){

	// Variáveis auxiliares do problema:
	constexpr auto gamma = g*dt/(rho*cp);
	constexpr auto r = kappa*dt/(rho*cp*dx*dx);

	// Malha:
	std::vector<std::vector<double>> T(Ny, std::vector<double>(Nx, T_init));
	// Arquivo de saida:
	std::fstream print_start {"2D_Heat.dat", std::ios::out|std::ios::trunc};
	std::fstream print_resume {"2D_Heat.dat", std::ios::out|std::ios::app};

	/*
	//Inicio da execução:

	//Manipulador do arquivo de saida:
	std::string arquivo {"Data_"};
	auto tempo = std::to_string(3.1415);
	std::fstream printer {arquivo + tempo + ".dat", std::ios::out|std::ios::trunc};
	std::fstream {"2D_Heat.dat", std::ios::out|std::ios::app};
	start_save_data(printer);
	*/
	std::cout << "\nExecution reached the end" << std::endl;
}
