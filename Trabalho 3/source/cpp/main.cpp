#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm> //std::fill
#include <fstream>   //std::fstream

// Protótipos de funções:
void solver_finite_difference(std::vector<std::vector<double>>& T);
void solver_variable_property(std::vector<std::vector<double>>& T);
void save_data_final(std::fstream& u_file, const std::vector<std::vector<double>>& T);
void fix_bounds(std::vector<std::vector<double>>& T, bool dirichlet = false);
void print_messages();
void try_static_assertions();
template <typename T> void print_array_2D(const std::vector<std::vector<T>> A);

// Variáveis do domínio e da simulaçaão:
constexpr double kappa {0.6};                           // coeficiente de difusividade térmica
constexpr double rho {600};                             // massa específica
constexpr double cp {1200};                             // calor específico
constexpr double h {15.0};                              // coeficiente de troca termica por convecao
constexpr double g {100000};                            // geração interna
constexpr double Lx {0.05};                             // tamanho em x
constexpr double Ly {0.05};                             // tamanho em y
constexpr int Nx {9};                                 // número de nós em x
constexpr int Ny {9};                                 // número de nós em y
constexpr auto dx = Lx/(Nx-1);
constexpr auto dy = Ly/(Ny-1);
constexpr double T_init {20.0};                         // temperatura inicial da placa
constexpr double T_fixed {100.0};                       // temperatura prescrita no contorno
constexpr double T_inf {20.0};                          // temperatura do fluido no interior
constexpr double ti {0.0};                              // tempo inicial da simulação
constexpr double tf {240};                              // tempo final da simulação
constexpr auto alfa = kappa/(rho*cp);
constexpr auto stab = (dx*dx)/(2*alfa);                 // máximo passo de tempo para estabilidade
constexpr auto dt = 0.95*stab;                          // passo de tempo (95 % do passo máximo)
constexpr int nsteps = 1 + static_cast<int>((tf-ti)/dt);// número de passos de tempo
constexpr auto lambda = g*dt/(rho*cp);
constexpr auto r = kappa*dt/(rho*cp*dx*dx);
constexpr auto gamma = (2*h*dx)/kappa;
auto p = static_cast<int>(Nx*0.4);                   // começo da área vazada: 30 % da extensão
auto p_w = static_cast<int>(Nx*0.6);                 // final da área vazada: 70 % da extensão


int main(int argc, char* argv[]){

	try_static_assertions();
	print_messages();

	// Malha:
	std::vector<std::vector<double>> T(Ny, std::vector<double>(Nx, T_init));
	// Arquivo de saida:
	std::fstream save_1 {"2D_Heat.txt", std::ios::out|std::ios::trunc};

	// Inicio da iteração temporal:
	for (int step = 0; step < 1; step++){
		fix_bounds(T);
		solver_variable_property(T);
		if (step%1000 == 0)
			std::cout << "Current progress: " << 100*step/nsteps << " %"<< std::endl;
	}
	save_data_final(save_1, T);

	std::cout << "\n------- Execution reached the end.------- \n" << std::endl;
}


void solver_finite_difference(std::vector<std::vector<double>>& T){

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
				T[i][j] = T[i][j] + r*(T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1] - 4*T[i][j]) + lambda;
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
	for (int i = p; i < p_w; i++){
		for (int j = p; j < p_w; j++){
			T[i][j] = T_inf;
		}
	}
}

double function_k1(double T){
	return 14.835 + 0.0174*T - 4.05*(10e-06) * (T*T);
}

void solver_variable_property(std::vector<std::vector<double>>& T){

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
				print_array_2D
				double coef = dt/(dx*dx*rho*cp);
				double T_ij = T[i][j];
				double Ti_prev = T[i-1][j];
				double Ti_next = T[i+1][j];
				double Tj_prev = T[i][j-1];
				double Tj_next = T[i][j+1];
				
				double ki_prev = (function_k1(T_ij) + function_k1(Ti_prev))/(2.0);
				double ki_next = (function_k1(T_ij) + function_k1(Ti_next))/(2.0);
				double kj_prev = (function_k1(T_ij) + function_k1(Tj_prev))/(2.0);
				double kj_next = (function_k1(T_ij) + function_k1(Tj_next))/(2.0);
				T[i][j] = T[i][j] + coef*(ki_next*T[i+1][j] + ki_prev*T[i-1][j] + kj_next*T[i][j+1] + kj_prev*T[i][j-1] - (ki_prev + ki_next + kj_prev + kj_next)*T[i][j]) + g*dt;
			}
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
	dt < stab ? std::cout << "STABILITY CRITERIA FULLFILED!" << std::endl : std::cout << "STABILITY CRITERIA VIOLATED!" << std::endl;
	std::cout << "Maximum time step allowed: " << stab << std::endl;
	std::cout << "Time step used: " << dt << std::endl; 
	std::cout << "dx = " << dx << std::endl;
	std::cout << "Coeficient r = " << r << std::endl;
	std::cout << "Total time steps to be evaluated: " << nsteps << std::endl;
}

void try_static_assertions(){

	static_assert(Nx != 0 && Ny != 0, "Numero de nos nao pode ser nulo");
	static_assert(Nx == Ny, "Numero de nos precisa em x e y ser igual");
	static_assert(Nx > 8 || Ny > 8, "Numero de nos precisa ser maior que 8");
}
template <typename T>
void print_array_2D(const std::vector<std::vector<T>> A){

	auto nrow = A.size();
	auto ncol = A[0].size();
	for(int i = 0; i < nrow; i++){
		for(int j = 0; j < ncol; j++){
			std::cout << std::setw(6) << A[i][j] << ' ';
		}
		std::cout << '\n';
	}
}