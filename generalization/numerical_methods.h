#pragma once

#include <vector>
#include <functional>


/*
du/dt + c du/dx = 0    c>0
*/

using vec1D = std::vector<double>;

struct WaveProblem
{
	double dx{0.1};
	double dt{0.1};
	double c{0.1};
	double tf{10.0};
	double length{1.0};
	double temp_left{1.0};
};


vec1D upwindFirstOrder(const WaveProblem& setup)
{
	auto nsteps{ setup.tf / setup.dt };                    // numero de passos para discretizar tf em intervalos de dt
	auto nx = static_cast<int>(setup.length/setup.dx) + 1; // numero de pontos para discretizar L em segmentos de tamanho dx
	vec1D u(nx, 0.0);

	u[0] = setup.temp_left;                                // valor prescrito a esquerda
	
	auto r{ setup.c*setup.dt/setup.dx};
	for (int n = 0; n < nsteps; ++n)                       // evolucao temporal
	{
		for (int j = 1; j < nx; ++j)                       // evolucao espacial
		{
			u[j] = r * (u[j-1] - u[j]) + u[j];
		}
	}
	return u;
}
vec1D laxMethod(const WaveProblem& setup)
{
	auto nsteps{ setup.tf / setup.dt };                      // numero de passos para discretizar tf em intervalos de dt
	auto nx = static_cast<int>(setup.length/setup.dx) + 1;   // numero de pontos para discretizar L em segmentos de tamanho dx
	vec1D u(nx, 0.0);

	u[0] = setup.temp_left;                                  // valor prescrito a esquerda
	u[nx-1] = setup.temp_left;                               // valor prescrito a direita
	
	auto r{ 0.5*setup.c * setup.dt / setup.dx };
	for (int n = 0; n < nsteps; ++n)                         // evolucao temporal
	{
		for (int j = 1; j < nx-1; ++j)                       // evolucao espacial
		{
			u[j] = r*(u[j-1] - u[j+1]) + 0.5*(u[j+1] - u[j-1]);
		}
	}
	
	return u;
}