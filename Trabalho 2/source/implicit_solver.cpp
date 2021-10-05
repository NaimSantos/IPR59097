#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <chrono>

#include <armadillo>

// Variáveis do problemna e do domínio da simulação:
constexpr double kappa{ 0.6 };                            // coeficiente de condutividade térmica
constexpr double h{ 15.0 };                               // coeficiente de troca de calor por convecção
constexpr double cp{ 1200.0 };                            // calor específico
constexpr double rho{ 600.0 };                            // massa específica
constexpr double g{ 100000.0 };                           // geração interna
constexpr double T_inf{ 20.0 };                           // temperatura do ambiente
constexpr double T0{ 20.0 };                              // temperatura inicial da placa
constexpr double L{ 0.3 };                                // comprimento da placa
constexpr int N{ 31 };                                    // número de nós da malha
constexpr double ti{ 0.0 };                               // tempo inicial da simulação
constexpr double tf{ 1000.0 };                            // tempo total da simulação
constexpr auto dx = L / (N - 1);                          // comprimento do intervalo
constexpr auto dt = 0.01;                                 // passo de tempo
constexpr int nsteps = static_cast<int>((tf - ti) / dt);  // número de passos de tempo
constexpr auto alpha = kappa / (rho * cp);
constexpr auto r1 = (alpha * dt) / (dx * dx);             // coeficiente r do método implícito
constexpr auto eta = 3.0 + ((2 * h * dx) / kappa);
constexpr auto sigma = (g * dt) / (rho * cp);

void implicit_diff(arma::mat& A, arma::vec& B, const double r);
void save_data(const arma::vec& Temperature, const arma::vec& Position, const int step);

int main(int argc, char** argv){
	std::cout << "\ndx = " << dx << ", dt = " << dt << std::endl;
	std::cout << "Calculating " << nsteps << " steps..." << std::endl;

	// Matrizes utilizadas:
	arma::mat A(N, N, arma::fill::zeros);
	arma::vec B(N, arma::fill::value(T0));
	arma::vec X = arma::linspace(0.0, L, N);

	implicit_diff(A, B, r1);

	save_data(B, X, static_cast<int>(tf));
	std::cout << "Execution reached the end" << std::endl;
	std::cin.get();
	return 0;
}

void implicit_diff(arma::mat& A, arma::vec& B, const double r){
	std::cout << "Numerical solver started..." << std::endl;
	// Preenchemos a matriz A (B já foi preenchido com os valores de T0):
	A(0,0) = -3.0;
	A(0,1) = 4.0;
	A(0,2) = -1.0;
	int k = 0;
	for (int i = 1; i < N-1; i++) {
		A(i, k) = -r;
		A(i, k+1) = 1 + 2 * r;
		A(i, k+2) = -r;
		k++;
	}
	A(N-1, N-3) = 1.0;
	A(N-1, N-2) = -4.0;
	A(N-1, N-1) = eta;

	// Os passos iterativos do método:
	for (int step = 1; step < nsteps; step++) {
		// Corrige B:
		B(0) = 0.0;
		for (int i = 1; i < N - 1; i++) {
			B(i) = B(i) + sigma;
		}
		B(N - 1) = (2 * h * dx * T0) / kappa;
		// Resolve o sistema:
		B = solve(A, B);
	}
	std::cout << "Numerical solver finished." << std::endl;
}
void save_data(const arma::vec& Temperature, const arma::vec& Position, const int step) {
	std::fstream printer{ "Temperatura_t_" + std::to_string(step) + ".txt", std::ios::out | std::ios::trunc };

	printer << "Pos Temperatura\n";
	for (int i = 0; i < N; i++) {
		printer << Position(i) << ' ' << Temperature(i) << '\n';
	}
}