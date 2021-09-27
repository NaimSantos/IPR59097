#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>

#include "utilities.h"

void fill_coef(std::vector<std::vector<double>>& A, const double gamma);
void fill_indepen(std::vector<double>& B, const double gamma);
void save_data(const std::vector<double>& V, const std::vector<double>& W, std::string& filename);
void save_analitic(const std::vector<double>& V1, const std::vector<double>& V2, std::string& filename);
double kelvin_to_celsius(const double TK);
// Variáveis do problema:
constexpr double kappa {59.0};                                           // coeficiente de condutividade térmica
constexpr double h {600};                                                // coeficiente de troca de calor por convecção
constexpr double D {4e-03};                                              // diâmetro da seção transversal da aleta
constexpr double T_amb{20.0};                                            // temperatura do ambiente
constexpr double T_0 {320};                                              // temperatura da aleta em x=0
constexpr double T_n {75};                                               // temperatura da aleta em x=L
constexpr double L {0.25};                                               // comprimento da aleta
constexpr int N {51};                                                    // número de nós em x

int main(int argc, char* argv[]){
	constexpr auto dx = L/(N-1);                                         // comprimento do intervalo em x
	constexpr auto A = NPI*D*D/4;                                        // área da seção transversal da aleta
	constexpr auto P = NPI*D;                                            // perímetro da aleta
	constexpr auto m2 = (h*P)/(kappa*A);
	constexpr auto gamma = dx*dx*m2;
	const auto m = std::sqrt(m2);
	std::cout << "dx = " << dx << std::endl;
	std::cout << "m2 = " << m2 << std::endl;
	std::cout << "m = " << std::sqrt(m2) << std::endl;
	std::cout << "m2*T_amb = " << m2*T_amb << std::endl;

	std::vector<std::vector<double>> G(N, std::vector<double>(N, 0.0));  // matriz NxN de coeficientes
	std::vector<double> T(N, 0.0);                                       // vetor onde serão armazenadas as temperaturas
	std::vector<double> B(N, 0.0);                                       // vetor b no sistema Ax=b
	std::vector<double> Pos(N, 0.0);                                     // vetor com as posições avaliadas
	std::vector<double> D(N, 0.0);                                       // vetor para analítico
	linspace(Pos, N, L, 0.0);

	// Corrigir os coeficientes e os termos independentes:
	fill_coef(G, gamma);
	fill_indepen(B, gamma);

	// Solução do sistema via Gauss-Siedel:
	GS_Solver<double>(G, B, T);

	std::string numeric_data {"numeric_data"};
	numeric_data += "_k_" + std::to_string(static_cast<int>(kappa)) + "_h_" + std::to_string(static_cast<int>(h)) + ".txt";
	save_data(T, Pos, numeric_data);



	// Solução analítica, para comparação posterior:
	for (int i = 0; i < N-1; i++){
		D[i] = T_amb + (T_0-T_amb) * (((T_n-T_amb)/(T_0-T_amb))*(std::sinh(m*Pos[i])) + std::sinh(m*(L - Pos[i]))) / (std::sinh(m*L));
	}
	std::string analitic_data {"data_analitic"};
	analitic_data += "_k_" + std::to_string(static_cast<int>(kappa)) + "_h_" + std::to_string(static_cast<int>(h)) + ".txt";
	save_analitic(Pos, D, analitic_data);


	std::cout << "Execution reached the end" << std::endl;
}

void fill_coef(std::vector<std::vector<double>>& A, const double gamma){
	int p = A.size();
	int q = A[0].size();
	if (p!=q || p!= N)
		return;

	A[0][0] = 1.0;
	for (int i = 1; i < (p-1); i++){
		A[i][i-1] = 1.0;
		A[i][i] = - (2 + gamma);
		A[i][i+1] = 1.0;
	}
	A[p-1][q-1] = 1.0;
}
void fill_indepen(std::vector<double>& B, const double gamma){
	auto q = B.size();
	if (q != N)
		return;

	B[0] = T_0;
	for (int i = 1; i < q-1; i++)
		B[i] = -gamma * T_amb;
	B[q-1] = T_n;
}
void save_data(const std::vector<double>& V, const std::vector<double>& W, std::string& filename){
	std::fstream saver {filename, std::ios::out|std::ios::trunc};
	saver << "Perfil de Temperatura\n";
	saver << "Posicao(m)\tTemperatura\n";

	const auto nelem = V.size();
	const auto nelem2 = W.size();
	if (nelem != nelem2)
		return;

	for (int i = 0; i < nelem; i++){
		saver << std::setw(8) << W[i] << "\t" << V[i] << "\n";
	}
}
void save_analitic(const std::vector<double>& V1, const std::vector<double>& V2, std::string& filename){
	std::fstream saver {filename, std::ios::out|std::ios::trunc};
	saver << "Solucao analitica\n";
	saver << "Posicao(m)\tTemperatura\n";

	const auto nelem = V1.size();
	const auto nelem2 = V2.size();
	if (nelem != nelem2)
		return;

	for (int i = 0; i < nelem; i++){
		saver << std::setw(8) << V1[i] << "\t" << V2[i] << "\n";
	}
}
double kelvin_to_celsius(const double TK){
	return TK - 273.15;
}
