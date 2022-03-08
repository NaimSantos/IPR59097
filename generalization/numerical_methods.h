#pragma once

#include <vector>
#include <functional>

/*
du/dt + c du/dx = 0    c>0
*/

using vec1D = std::vector<double>;
using vec2D = std::vector<std::vector<double>>;

struct WaveProblem
{
	double length{1.0};     // tamanho do dominio
	double tf{10.0};        // tempo final de simulacao
	double dx{0.1};
	double dt{0.1};
	double c{0.1};          // velocidade de propagacao

	double temp_left{1.0};
};


vec1D upwindFirstOrder(const WaveProblem& input)
{
	double nsteps{input.tf/input.dt};                        // numero de passos para discretizar tf em intervalos de dt
	int nx = static_cast<int>(input.length/input.dx) + 1; // numero de pontos para discretizar L em segmentos de tamanho dx
	vec1D u(nx, 0.0);

	// Upwind requer que a velocidade nao seja nula:
	if (c = 0.0)
		return u;

	auto CFL{input.c*input.dt/input.dx};
	// Estabilidade: CFL <= 1.0
	if (CFL > 1.0)
		return u;

	// Tratamento de cada caso do sinal da velocidade:
	if (c > 0.0)
	{
		for (int n = 0; n < nsteps; ++n)
		{
			for (int j = 1; j < nx; ++j)
			{
				u[j] = CFL*(u[j-1] - u[j]) + u[j];
			}
		}
	}
	else
	{
		for (int n = 0; n < nsteps; ++n)
		{
			for (int j = 1; j < nx; ++j)
			{
				u[j] = CFL*(u[j] - u[j+1]) + u[j];
			}
		}
	}

	return u;
}
vec1D laxMethod(const WaveProblem& input)
{
	auto nsteps{ input.tf / input.dt };                      // numero de passos para discretizar tf em intervalos de dt
	auto nx = static_cast<int>(input.length/input.dx) + 1;   // numero de pontos para discretizar L em segmentos de tamanho dx
	vec1D u(nx, 0.0);

	u[0] = input.temp_left;                                  // valor prescrito a esquerda
	u[nx-1] = input.temp_left;                               // valor prescrito a direita
	
	auto r{ 0.5*input.c * input.dt / input.dx };
	for (int n = 0; n < nsteps; ++n)                         // evolucao temporal
	{
		for (int j = 1; j < nx-1; ++j)                       // evolucao espacial
		{
			u[j] = r*(u[j-1] - u[j+1]) + 0.5*(u[j+1] - u[j-1]);
		}
	}
	
	return u;
}