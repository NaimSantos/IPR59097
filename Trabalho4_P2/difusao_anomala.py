# APÊNDICE B - DIFUSÃO ANÔMALA

import numpy as np
import matplotlib.pyplot as plt
import math

# Variáveis do problema e do domínio:
L = 1.0                  # comprimento total do domínio
N  = 21                  # número de nós da malha
ti = 0.0                 # tempo inicial da simulação
tf = 10.0                # tempo final da simulação
dx = L/(N-1)             # comprimento do intervalo
dt = 0.001               # passo de tempo
nsteps = int((tf-ti)/dt) # número de passos de tempo

# parâmetros do artigo de Silva et al.:
K2 = 0.001
K4 = 0.00001
beta = 0.2
# parâmetros do artigo do Sulivan: beta = 0.5; K2 = 0.1; K4 = 0.00001

kappa1 = beta*K2*dt/(dx**2)
kappa2 = - ((beta*(1-beta))*K4*dt)/(dx**2)
a1 = -kappa2
a2 = ((4*kappa2) - kappa1)
a3 = (1+(2*kappa1) - (6*kappa2))

# Malhas:
X = np.linspace(0.0, L, N)    # posições, para plotar
Analytic = np.zeros(N)        # Solução analítica
A =  np.zeros((N, N))         # Matriz dos coeficientes
B = np.zeros(N)               # Matriz dos termos independentes

def plot_perfil_single(x, y):
    plt.plot(x, y, 'r', linewidth=2)
    plt.xlabel("Posição (m)", fontsize = 11)
    plt.ylabel("Concentração", fontsize = 11)
    plt.grid(True, 'major', 'both')
    plt.savefig('Concentracao_no_Espaco.png')
    plt.show()

def plot_compare(x, y1, y2):
    title = "Solução em t = {:.2f} s".format(tf)
    plt.title(title)
    plt.plot(x, y1, 'b', label='Analítico', linewidth=2)
    plt.plot(x, y2, 'r', label='Numérico', linestyle='dashed', linewidth=2)
    plt.title(title)
    plt.xlabel("Posição (m)", fontsize = 11)
    plt.ylabel("Concentração", fontsize = 11)
    plt.grid(True, 'major', 'both')
    plt.legend(loc='best', fontsize=8)
    plt.savefig('Comparacao_Analitico_Numerico.png')
    plt.show()


def initial_condition(x):
    return math.sin(math.pi*x/L)

def solucao_analitica(x, t):
    return math.exp(-beta*K2*(math.pi**2)*t)*math.exp(-beta*(1-beta)*K4*(math.pi**4)*t)*math.sin(math.pi*x/L)

def solve_explicitly():
    global A
    global B

    # Preenchemos A:
    # Primeira linha:
    A[0, 0] = 1.0
    # Segunda linha:
    A[1, 0] = 2.0; A[1, 1] = -5.0; A[1, 2] = 4.0; A[1, 3] = -1.0;
    # Da terceira (índice 2) até a antepenúltima (índice N-3) linha:
    for i in range(2, N-2):
        A[i, i-2] = a1
        A[i, i-1] = a2
        A[i, i]   = a3
        A[i, i+1] = a2
        A[i, i+2] = a1
    # Penúltima linha:
    A[N-2, N-4] = -1.0;  A[N-2, N-3] = 4.0; A[N-2, N-2] = -5.0; A[N-2, N-1] = 2.0;
    # Ultima linha:
    A[N-1, N-1] = 1.0

    # O vetor de termos independentes é, inicialmente, os P no instante 0, então preenchemos B com a condição inicial
    # Preenchemos da terceira (índice 2) até a antepenúltima (índice N-3) célula:
    for i in range (2, N-2):
        B[i] = initial_condition(i*dx)

    # Iteração no tempo:
    for n in range(0, nsteps):
        # Corrigimos os Bs:
        B[0] = 0.0
        B[1] = 0.0
        B[N-2] = 0.0
        B[N-1] = 0.0
        # Obtemos o novo B como a solução do sistema linear:
        B = np.linalg.solve(A, B)


# Solução analítica:
for i in range (0, N):
    Analytic[i] = solucao_analitica(i*dx, tf)


# Solução numérica:
solve_explicitly()
plot_compare(X, Analytic, B)

