#        APÊNDICE A - SOLUÇÃO IMPLÍCITA

import numpy as np
import matplotlib.pyplot as plt
import math

# Variáveis do problema e do domínio:
L = 0.03                 # comprimento total da placa
N  = 51                  # número de nós da malha
ti = 0.0                 # tempo inicial da simulação
tf = 12000.0             # tempo final da simulação
dx = L/(N-1)             # comprimento do intervalo
dt = 0.1                 # passo de tempo
nsteps = int((tf-ti)/dt) # número de passos de tempo
kappa = 0.6
rho = 600
cp = 1200
h = 15.0
g = 100000
T0 = 20.0
TL = 20.0
r1 = (kappa*dt)/(rho*cp*dx*dx) 
eta = 3.0 + ((2*h*dx)/kappa)
mu = (dt)/(rho*cp)

print("Refinamento da malha espacial: ", dx)
print("Criterio de estabiliade: ", r1)
print("Passo de tempo: ", dt)
print("Numero de passos de tempo: ", nsteps)
T1 = np.zeros((nsteps, N))           # para armazenar os resultados em todos os tempos
X = np.linspace(0.0, L, N)           # posições, para plotar
tempos = np.linspace(ti, tf, nsteps) # tempos, para plotar

def function_g(x):
    return math.exp(x)
    
def plot_perfil_single(x, y):
    plt.plot(x, y, 'r', linestyle='dashed', linewidth=2)
    plt.xlabel("Posição (m)", fontsize = 11)
    plt.ylabel("Temperatura (° C)", fontsize = 11)
    #plt.legend(loc='upper center', fontsize=9)
    plt.grid(True, 'major', 'both')
    plt.savefig('Grafico_Perfil.png')
    plt.show()

def plot_perfil_tri(x, y1, y2, y3):
    plt.plot(x, y1, 'c', label= 't = 1 s',  linewidth=2,)
    plt.plot(x, y2, 'b', label= 't = 1000 s', linewidth=2, )
    plt.plot(x, y3, 'r', label= 't = 120000 s', linewidth=2,)
    plt.xlabel("Posição (m)", fontsize = 11)
    plt.ylabel("Temperatura (° C)", fontsize = 11)
    plt.legend(loc='best', fontsize=8)
    plt.grid(True, 'major', 'both')
    plt.savefig('Grafico_Perfil_3_Tempos.png')
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

def implictsolver(A, B, T):
    T[0] = B.reshape(1, N)
    t = 1
    while t < nsteps:
        # Correção de B:
        B[0][0] = 0.0
        B[N-1][0] = (2*h*dx*T0)/kappa
        i = 1
        while i < N-1 :
            B[i][0] = B[i][0] + mu*g
            #B[i][0] = B[i][0] + mu*function_g(i*dx)*100000
            i = i + 1
        # Atualiza B com as novas temperaturas:
        B = np.linalg.solve(A, B)
        T[t] = B.reshape(1, N)
        t = t + 1

def solveimplicitly(r, T):
    # Preenchimento da matriz de termos independentes:
    B = np.full((N, 1), T0)
    # Preenchimento da matriz de coeficientes:
    A = np.zeros((N,N))
    A[0][0] = -3.0
    A[0][1] = 4.0
    A[0][2] = -1.0
    A[N-1][N-3] = 1.0
    A[N-1][N-2] = -4.0
    A[N-1][N-1] = eta
    i = 1
    j = 0
    while i < N-1 :
        A[i][j] = - r
        A[i][j+1] = 1 + 2*r
        A[i][j+2] = -r
        j = j+1
        i = i+1
    implictsolver(A, B, T)

# O método é chamado aqui, B será alterado:
solveimplicitly(r1, T1)

# Apenas para encontrar os indices de x=0.25L, 0.5L e 0.75L
index = np.linspace(0, N, 5)
index_x1 = (int)(index[1] - 1)
index_x2 = (int)(index[2] - 1)
index_x3 = (int)(index[3] - 1)
Y1=T1[:, index_x1]
Y2=T1[:, index_x2]
Y3=T1[:, index_x3]

# Apenas para encontrar o indices dos tempos requeridos:
i1 = (int)(1/dt)
i120 = (int)(1000/dt)

# Plota os 3 gráficos de resultados:
plot_evolution(tempos, Y1, Y2, Y3)
plot_perfil_tri(X, T1[i1-1], T1[i120-1], T1[nsteps-1])
plot_perfil_single(X, T1[nsteps-1])