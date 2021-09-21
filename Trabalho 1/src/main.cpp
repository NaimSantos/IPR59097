#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>

#include "utilities.h"

// Protótipos de funções:
void fill_coef(std::vector<std::vector<double>>& A, const double gamma);
void save_data(const std::vector<double>& V, const std::vector<double>& W);
double kelvin_to_celsius(const double TC);
// Variáveis do problema:
constexpr double kappa {59.0};                        // coeficiente de condutividade térmica
constexpr double h {380};                             // coeficiente de troca de calor por convecção
constexpr double D {4e-03};                           // diâmetro da seção transversal da aleta
constexpr double T_amb{20.0};                         // temperatura do ambiente
constexpr double T_0 {320};                         // temperatura da aleta em x=0
constexpr double T_n {75};                          // temperatura da aleta em x=L
constexpr double L {0.25};                            // comprimento da aleta
constexpr int N {50};                                // número de nós em x

int main(int argc, char* argv[]){
	constexpr auto dx = L/(N-1);                                             // comprimento do intervalo em x
	constexpr auto A = NPI*D*D/4;                                            // área da seção transversal da aleta
	constexpr auto P = NPI*D;                                                // perímetro da aleta
	constexpr auto m2 = (h*P)/(kappa*A);
	constexpr auto gamma = dx*dx*m2;
	const auto m = std::sqrt(m2);
	std::cout << "dx = " << dx << std::endl;
	std::cout << "m2 = " << m2 << std::endl;
	std::cout << "m = " << std::sqrt(m2) << std::endl;
	std::cout << "m2*T_amb = " << m2*T_amb << std::endl;

	std::vector<std::vector<double>> G(N, std::vector<double>(N, 0.0));      // matriz NxN de coeficientes
	std::vector<double> T(N, 0.0);                                           // vetor onde serão armazenadas as temperaturas
	std::vector<double> B(N, 0.0);                                           // vetor b no sistema Ax=b
	std::vector<double> Pos(N, 0.0);                                         // vetor com as posições avaliadas
	linspace(Pos, N, L, 0.0);

	fill_coef(G, gamma);

	B[0] = T_0;
	for (int i = 1; i < N-1; i++)
		B[i] = -gamma * T_amb;
	B[N-1] = T_n;

	GS_Solver<double>(G, B, T);
	save_data(T, Pos);

	std::cout << "Execution reached the end" << std::endl;
}

void fill_coef(std::vector<std::vector<double>>& A, const double gamma){
	int m = A.size();
	int n = A[0].size();

	A[0][0] = 1.0;
	for (int i = 1; i < (m-1); i++){
		A[i][i-1] = 1.0;
		A[i][i] = - (2 + gamma);
		A[i][i+1] = 1.0;
	}
	A[m-1][n-1] = 1.0;

}
void save_data(const std::vector<double>& V, const std::vector<double>& W){
	std::fstream saver {"temperature_out1.dat", std::ios::out|std::ios::trunc};
	saver << "Perfil de Temperatura\n";
	saver << "Posicao(m)\tTemperatura\n";

	const auto N = V.size();
	for (int i = 0; i < N; i++){
		saver << std::setw(8) << W[i] << "\t" << V[i] << "\n";
	}
}

double kelvin_to_celsius(const double TC){
	return TC - 273.15;
}