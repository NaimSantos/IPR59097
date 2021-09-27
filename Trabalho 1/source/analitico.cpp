#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include "utilities.h"

int main (int argc, char* argv[]){
	constexpr double kappa {59.0};                        // coeficiente de condutividade térmica
	constexpr double h {380};                             // coeficiente de troca de calor por convecção
	constexpr double D {4e-03};                           // diâmetro da seção transversal da aleta
	constexpr double T_amb{20.0};                         // temperatura do ambiente
	constexpr double T_0 {320};                           // temperatura da aleta em x=0
	constexpr double T_n {75};                            // temperatura da aleta em x=L
	constexpr double L {0.25};                            // comprimento da aleta
	constexpr int N {51};                                 // número de nós
	constexpr auto dx = L/(N-1);                          // comprimento do intervalo em x
	constexpr auto A = NPI*D*D/4;                         // área da seção transversal da aleta
	constexpr auto P = NPI*D;                             // perímetro da aleta
	constexpr auto m2 = (h*P)/(kappa*A);
	constexpr auto gamma = dx*dx*m2;
	const auto m = std::sqrt(m2);
	std::cout << "dx = " << dx << std::endl;
	std::cout << "m2 = " << m2 << std::endl;
	std::cout << "m = " << std::sqrt(m2) << std::endl;
	std::cout << "m2*T_amb = " << m2*T_amb << std::endl;

	std::vector<double> X(N, 0.0);                       // vetor com as posições avaliadas
	linspace(X, N, L, 0.0);
	std::vector<double> B(N, 0.0);
	for (int i = 0; i < N-1; i++){
		B[i] = T_amb + (T_0-T_amb) * (((T_n-T_amb)/(T_0-T_amb))*(std::sinh(m*X[i])) + std::sinh(m*(L - X[i]))) / (std::sinh(m*L));
	}
	std::fstream saver {"analitico.txt", std::ios::out|std::ios::trunc};
	saver << "Perfil de Temperatura\n";
	saver << "Posicao(m)\tTemperatura\n";
	for (int i = 0; i < N; i++){
		saver << std::setw(8) << X[i] << "\t" << B[i] << "\n";
	}
}