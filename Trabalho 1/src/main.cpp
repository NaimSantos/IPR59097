#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>

#include "utilities.h"

// Protótipos de funções:
void fill_coef(std::vector<std::vector<double>>& A, const double gamma);
void fill_b(std::vector<double>& A, const double gamma);
void save_data(const std::vector<double>& V, const std::vector<double>& W);

// Variáveis do problema:
constexpr double kappa {59.0};                        // coeficiente de condutividade térmica
constexpr double h {380};                             // coeficiente de troca de calor por convecção
constexpr double D {4e-03};                           // diâmetro da seção trasnversal da aleta
constexpr double T_amb{20.0};                         // temperatura do ambiente
constexpr double T_0 {320.0};                         // temperatura da aleta em x=0
constexpr double T_n {75.0};                          // temperatura da aleta em x=L
constexpr double L {0.25};                            // comprimento da aleta
constexpr int N {10};                                 // número de nós em x  


int main(int argc, char* argv[]){

	constexpr auto A = NPI*D*D/4;                                            // área da seção transversal da aleta
	constexpr auto P = NPI*D;                                                // perímetro da aleta
	constexpr auto m2 = (h*P)/(kappa*A);
	constexpr auto dx = L/(N-1);                                             // comprimento do intervalo em x
	constexpr auto gamma = dx*m2;

	std::vector<std::vector<double>> G(N, std::vector<double>(N, 0.0));      // matriz de coeficientes NxN, inicializada em 0.0
	std::vector<double> X(N, 0.0);                                           // vetor onde serão armazenadas as temperaturas
	std::vector<double> B(N, 0.0);                                           // vetor b no sistema Ax=b
	std::vector<double> Pos(N, 0.0);                                         // vetor com as posições avaliadas
	linspace(Pos, N, L, 0.0);

	fill_coef(G, gamma);
	fill_b(B, gamma);
	GS_Solver(G, B, X);
	save_data(X, Pos);

	std::cout << "\nExecution reached the end" << std::endl;
}
void fill_coef(std::vector<std::vector<double>>& A, const double gamma){
	int m = A.size();
	int n = A[0].size();
	// Corrige da primeira até a penúltima linha
	for (int i = 1; i < m-1; i++){
		A[i][i-1] = 1.0;
		A[i][i]   = - (2 + gamma);
		A[i][i+1] = 1.0;
	}
	A[0][0] = 1.0;
	A[m][n] = 1.0;
}
void fill_b(std::vector<double>& A, const double gamma){
	auto m = A.size();
	for (auto& x : A){
		x = -gamma * T_amb;
	}
	A[0] = T_0;
	A[m-1] = T_n;
}

void save_data(const std::vector<double>& V, const std::vector<double>& W){
	// Arquivo de saida:
	std::fstream saver {"temperature_out.txt", std::ios::out|std::ios::trunc};
	
	saver << "Perfil de Temperatura\n";
	saver << "Posicao(m)\tTemperatura\n";

	const auto N = V.size();
	for (int i = 0; i < N; i++){
		saver << std::setprecision(2) << std::setw(8) << V[i] << "\t" << W[i] << "\n";
	}
}
