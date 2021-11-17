#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm> //std::fill
#include <fstream>   //std::fstream

void fix_bounds(std::vector<std::vector<double>>& T, int cond_type);
double g_ij(const double T);
double k_at_T(const double T);
void solver_finite_difference(std::vector<std::vector<double>>& T);
void save_data_final(std::fstream& u_file, const std::vector<std::vector<double>>& T);
void print_messages();
void try_static_assertions();

// Variáveis do domínio e da simulação:
constexpr double kappa {113};
constexpr double rho {7140};
constexpr double cp {3820};
constexpr double h {15.0};
constexpr double g {100000};
constexpr double Lx {0.05};
constexpr double Ly {0.05};
constexpr int Nx {201};
constexpr int Ny {201};
constexpr auto dx = Lx/(Nx-1);
constexpr auto dy = Ly/(Ny-1);
constexpr double T_init {20.0};
constexpr double T_fixed {100.0};
constexpr double T_fluid {20.0};
constexpr double ti {0.0};
constexpr double tf {120};  
constexpr auto alfa = kappa/(rho*cp);
constexpr auto stab = (dx*dx)/(2*alfa);
constexpr auto dt = 0.95*stab;
constexpr int nsteps = 1 + static_cast<int>((tf-ti)/dt);
constexpr auto lambda = dt/(rho*cp);
constexpr auto r = kappa*dt/(rho*cp*dx*dx);
constexpr auto gamma = (2*h*dx)/kappa;
auto p = static_cast<int>(Nx*0.3);
auto p_w = static_cast<int>(Nx*0.7);

int main(int argc, char* argv[]){
	try_static_assertions();
	print_messages();

	std::vector<std::vector<double>> T(Ny, std::vector<double>(Nx, T_init));
	std::fstream save_1 {"2D_Heat.txt", std::ios::out|std::ios::trunc};

	for (int step = 0; step < nsteps; step++){
		fix_bounds(T, 2);
		solver_finite_difference(T);
		if (step%1000 == 0)
			std::cout << "Current time step: " << step << std::endl;
	}
	save_data_final(save_1, T);
	std::cout << "\n------- Execution reached the end.------- ";
}

void fix_bounds(std::vector<std::vector<double>>& T, int cond_type){
	for (int i = p; i < p_w; i++){
		for (int j = p; j < p_w; j++){
			T[i][j] = T_fluid;
		}
	}
	for (int i = 0; i < Nx; i++){
		T[0][i] = T_fixed;
		T[i][Nx-1] = T_fixed;
	}
	switch (cond_type) {
		case 1 :
			for (int i = 0; i < Nx; i++){
				T[i][0] = T_fixed;
				T[Nx-1][i] = T_fixed;
			}
			break;
		case 2 :
			for (int i = 0; i < Nx; i++){
				T[i][0] = T[i][1];
				T[Nx-1][i] = T[Nx-2][i];
			}
			break;
		case 3 :
			break;
	}
}

double g_ij(const double T){
	return 8000*std::exp(-0.1*T);
}

double k_at_T(const double T){
	return 111.2 + 0.004*T + 0.000097*T*T;
}

double k_at_n(const double Ta, const double Tb){
	return 1.0/(0.5 *(1.0/Ta + 1.0/Tb));
}

void solver_finite_difference(std::vector<std::vector<double>>& T){


	for (int i = 1; i < Nx-1; i++){
		for (int j = 1; j < Ny-1; j++){
			if ((i>p && i<p_w) && (j>p && j<p_w)){
				continue;
			}
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
			else{
				T[i][j] = T[i][j] + r*(T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1] - 4*T[i][j]) + g_ij(T[i][j])*lambda;
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
	std::cout << "\nTotal time steps to be evaluated: " << nsteps << std::endl;
}

void try_static_assertions(){
	static_assert(Nx != 0 && Ny != 0, "Numero de nos nao pode ser nulo");
	static_assert(Nx == Ny, "Numero de nos precisa em x e y ser igual");
	static_assert(Nx > 10 || Ny > 10, "Numero de nos precisa ser maior que 10");
}
