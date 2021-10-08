import numpy as np
import matplotlib.pyplot as plt

# Variáveis do problema e do domínio:
L = 0.03                 # comprimento total da placa
N  = 25                  # número de nós da malha
ti = 0.0                 # tempo inicial da simulação
tf = 100.0               # tempo final da simulação
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
r = (kappa*dt)/(rho*cp*dx*dx) 
gamma = 3.0 + ((2*h*dx)/kappa)
mu = (dt)/(rho*cp)

A = np.zeros((N,N))                  # Matriz de coeficientes
B = np.full((N, 1), T0)              # Termos independente
X = np.linspace(0.0, L, N)           # Posições
tim = np.linspace(ti, tf, nsteps)    # Tempos, para plotar

def show_parameters():
    print('\nDomínio: %f (m)' %L)
    print('Refimamento: %.2f (m)' %dx)
    print("Critério de estabilidade: ", r)
    print('Tempo total de simulação: %0.1f (s)' %tf)
    print('Passo de tempo: %.1f (s)' %dt)
    print("Número de passos de tempo: ", nsteps)

def plot_perfil_single(x, y):
    plt.plot(x, y, 'r', linestyle='dashed', linewidth=2)
    plt.xlabel("Posição (m)", fontsize = 11)
    plt.ylabel("Temperatura (° C)", fontsize = 11)
    plt.grid(True, 'major', 'both')
    plt.savefig('Perfil_de_Temperatura.png')
    plt.show()

def solve_implicitly(nsteps):
    global A
    global B
    # Preenche A:
    A[0][0] = -3.0
    A[0][1] = 4.0
    A[0][2] = -1.0
    A[N-1][N-3] = 1.0
    A[N-1][N-2] = -4.0
    A[N-1][N-1] = gamma
    i, j = 1, 0
    while (i < N-1) :
        A[i][j] = - r
        A[i][j+1] = 1 + 2*r
        A[i][j+2] = -r
        j = j+1
        i = i+1
    # Itera de 1 até nsteps:
    t = 1
    while (t <= nsteps):
        # Correção dos B:
        B[0][0] = 0.0
        k = 1
        while k < (N-1) :
            B[k][0] = B[k][0] + mu*g
            k = k + 1
        B[N-1][0] = (2*h*dx*T0)/kappa
        # Atualiza B com as novas temperaturas:
        B = np.linalg.solve(A, B)
        t = t + 1
    print(t)

show_parameters()
solve_implicitly(nsteps)
plot_perfil_single(X, B)