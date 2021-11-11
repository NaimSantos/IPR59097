import numpy as np
import matplotlib.pyplot as plt
import math

# Variáveis do problema e do domínio:
L = 0.03                 # comprimento total da placa
N  = 51                  # número de nós da malha
ti = 0.0                 # tempo inicial da simulação
tf = 12000               # tempo final da simulação
dx = L/(N-1)             # comprimento do intervalo
dt = 0.2                 # passo de tempo
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
tempos = np.linspace(ti, tf, nsteps) # tempos, para plotar
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
def plot_evolution(x, y1, y2, y3):
    plt.plot(x, y1, 'r', label= 'x = 0.25L', linewidth=2)
    plt.plot(x, y2, 'b', label= 'x = 0.5L', linewidth=2)
    plt.plot(x, y3, 'g', label= 'x = 0.75L', linewidth=2)
    plt.xlabel("Tempo (s)", fontsize = 11)
    plt.ylabel("Temperatura (° C)", fontsize = 11)
    plt.legend(loc='lower right', fontsize=9)
    plt.grid(True, 'major', 'both')
    plt.savefig('Graf_Evolucao_Temporal.png')
    plt.show()

def solve_explicitly(T, nsteps):
    for i in range (1, nsteps):
        if (i%1000 == 0):
            print("Time step: ", i)
        # contorno esquerdo:
        T[i, 0] = T[i-1, 0] + 2*r*(T[i-1, 1] - T[i-1, 0]) + gamma*g
        # nós internos:
        for j in range (1, N-1):
            T[i, j] = T[i-1, j] + r*(T[i-1, j-1] - 2*T[i-1, j] + T[i-1, j+1]) + gamma*g
        # contorno direito:
        T[i, N-1] = T[i-1, N-1] + r*(2*T[i-1, N-2] - (2+mu)*T[i-1, N-1] + mu*T0) + gamma*g


show_parameters()
solve_explicitly(T, nsteps)
maxT = np.amax(T[nsteps-1])

# Apenas para encontrar os indices de x=0.25L, 0.5L e 0.75L
index = np.linspace(0, N, 5)
index_x1 = (int)(index[1] - 1)
index_x2 = (int)(index[2] - 1)
index_x3 = (int)(index[3] - 1)
Y1=T[:, index_x1]
Y2=T[:, index_x2]
Y3=T[:, index_x3]

# Plota o perfil final:
plot_perfil_single(X, T[nsteps-1, :])
# Plota os 3 gráficos de resultados:
plot_evolution(tempos, Y1, Y2, Y3)