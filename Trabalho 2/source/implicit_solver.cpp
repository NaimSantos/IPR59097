#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <chrono>

#include <armadillo>

// Variáveis do problemna e do domínio da simulação: :
constexpr double kappa{ 0.6 };                            // coeficiente de condutividade térmica
constexpr double h{ 15.0 };                               // coeficiente de troca de calor por convecção
constexpr double cp{ 1200.0 };                            // calor específico
constexpr double rho{ 600.0 };                            // massa específica
constexpr double g{ 100000.0 };                           // geração interna
constexpr double T_inf{ 20.0 };                           // temperatura do ambiente
constexpr double T0{ 20.0 };                              // temperatura inicial da placa
constexpr double L{ 0.3 };                                // comprimento da placa
constexpr int N{ 11 };                                    // número de nós da malha
constexpr double ti{ 0.0 };                               // tempo inicial da simulação
constexpr double tf{ 1000.0 };                            // tempo total da simulação
constexpr auto dx = L / (N - 1);                          // comprimento do intervalo
constexpr auto dt = 0.01;                                 // passo de tempo
constexpr int nsteps = static_cast<int>((tf - ti) / dt);  // número de passos de tempo
constexpr auto alpha = kappa / (rho * cp);
constexpr auto r1 = (alpha * dt) / (dx * dx);             // coeficiente r do método implícito
constexpr auto eta = 3.0 + ((2 * h * dx) / kappa);
constexpr auto sigma = (g * dt) / (rho * cp);

void implicit_diff(arma::mat& A, arma::mat& B, const double r);

int main(int argc, char** argv){
	std::cout << "\ndx = " << dx << ", dt = " << dt << std::endl;
	std::cout << "Calculating " << nsteps << " steps..." << std::endl;

	// Matrizes utilizadas:
	arma::mat A(N, N, arma::fill::zeros);
	arma::vec B(N, arma::fill::value(T0));

	arma::mat X = arma::linspace(0.0, L);
	X.print("X:");

	B.print("B:");
	implicit_diff(A, B, r1);
	B.print("B:");

	std::cout << "Execution reached the end" << std::endl;
	std::cin.get();
	return 0;
}

void implicit_diff(arma::mat& A, arma::mat& B, const double r){

	std::fstream printer{ "Temperatura_Implicit_Diff.txt", std::ios::out | std::ios::trunc };
	printer << "Perfil de Temperatura via Solver Implicito.\n";
	printer << "i t x T\n";

	int step = 0;
	printer << step << ' ' << ti << ' ';
	for (double x = 0.0; x <= L; x = x+dx)
		printer << x << ' ';
	printer << T0 << '\n';

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
	for (step = 1; step < nsteps; step++) {
		// Corrige B:
		B(0) = 0.0;
		for (int i = 1; i < N - 1; i++) {
			B(i) = B(i) + sigma;
		}
		B(N - 1) = (2 * h * dx * T0) / kappa;
		// Resolve o sistema:
		B = solve(A, B);
		printer << ' ' << step * dt << ' ';
		for (int i = 0; i < N; i++) {
			printer << B[i] << ' ';
		}
		printer << '\n';
		
	}
	printer.close();
}
