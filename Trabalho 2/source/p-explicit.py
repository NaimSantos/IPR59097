import numpy as np
import matplotlib.pyplot as plt
import math

# Variáveis do problema e do domínio:
L = 0.03                 # comprimento total da placa
N  = 31                  # número de nós da malha
ti = 0.0                 # tempo inicial da simulação
tf = 1500               # tempo final da simulação
dx = L/(N-1)             # comprimento do intervalo
dt = 0.5                 # passo de tempo
nsteps = int((tf-ti)/dt) # número de passos de tempo
kappa = 0.6
rho = 600
cp = 1200
h = 15.0
g = 100000
T0 = 20.0
TL = 20.0
alph = kappa/(rho*cp)
r = (kappa*dt)/(rho*cp*dx*dx)
gamma = alph*dt/kappa
mu = (2*h*dx)/kappa


X = np.linspace(0.0, L, N)           # Posições
tim = np.linspace(ti, tf, nsteps)    # Tempos, para plotar
T = np.full((nsteps, N), T0)         # Todos os tempos
def show_parameters():
    print('\nDomínio: %f (m)' %L)
    print('Refimamento: %.2f (m)' %dx)
    print("Critério de estabilidade: ", r)
    print('Tempo total de simulação: %0.1f (s)' %tf)
    print('Passo de tempo: %.1f (s)' %dt)
    print("Número de passos de tempo: ", nsteps)

def plot_perfil_single(x, y):
    plt.plot(x, y, 'r', linewidth=2)
    plt.xlabel("Posição (m)", fontsize = 11)
    plt.ylabel("Temperatura (° C)", fontsize = 11)
    plt.grid(True, 'major', 'both')
    plt.savefig('Perfil_de_Temperatura.png')
    plt.show()

def solve_explicitly(T, nsteps):
    for i in range (1, nsteps):
        print("Time Step =", i)
        T[i, 0] = T[i-1, 0] + 2*r*(T[i-1, 1] - T[i-1, 0]) + gamma*g
        for j in range (1, N-1):
            T[i, j] = T[i-1, j] + r*(T[i-1, j-1] - 2*T[i-1, j] + T[i-1, j+1]) + gamma*g
        T[i, N-1] = T[i-1, N-1] + r*(2*T[i-1, N-2] - (2+mu)*T[i-1, N-1] + mu*T0) + gamma*g


show_parameters()
solve_explicitly(T, nsteps)
plot_perfil_single(X, T[nsteps-1, :])