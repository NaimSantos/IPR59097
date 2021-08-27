#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <memory>    //std::unique_ptr
#include <algorithm> //std::fill
#include <fstream>   //std::fstream

#include "utilities.h"
#include "heat_transfer_2D.h"

//Variáveis do domímio e da simulação:
constexpr double kappa {0.6};                          // coeficiente de difusividade térmica
constexpr double T_init {20.0};							   // temperatura inicial da placa	
constexpr double L {1.0};                              // tamanho em x
constexpr double H {1.0};                              // tamanho em y
constexpr int Nx {10};                                 // número de nós em x
constexpr int Ny {10};                                 // número de nós em y
constexpr double ti {0.0};                             // tempo inicial da simulação
constexpr double tf {100.0};                          // tempo final da simulação
constexpr auto dx  = L/(Nx-1);                         // comprimento da célula
constexpr auto dy  = L/(Ny-1);                         // altura da célula
constexpr int nsteps {4096};                           // número de passos de tempo
constexpr auto dt = (tf-ti)/nsteps;                    // tamanho do passo de tempo

int main(int argc, char* argv[]){
	
	//Malha:
	std::vector<std::vector<double>> T(Ny, std::vector<double>(Nx, T_init));
	std::vector<double> T0 (Ny, 6.6);
	
	//Manipulador do arquivo de saida:
	std::string arquivo {"Data_"};
	auto tempo = std::to_string(3.1415);
	std::fstream printer {arquivo + tempo + ".dat", std::ios::out|std::ios::trunc};
	//std::fstream {"2D_Heat.dat", std::ios::out|std::ios::app};
	start_save_data(printer);
	std::cout << "\nExecution reached the end" << std::endl;
}