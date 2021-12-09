# APÊNDICE B - DIFUSÃO ANÔMALA

import numpy as np
import matplotlib.pyplot as plt
import math

# Variáveis do problema e do domínio:
L = 1.0                  # comprimento total do domínio
N  = 101                 # número de nós da malha
dx = L/(N-1)             # comprimento do intervalo
dt = 0.001                # passo de tempo
# Tempos de interesse:
tf1 = 0.5
tf2 = 4.0
tf3 = 10.0
# Parâmetros do modelo (artigo de Silva et al., 2014):
#beta = 0.2
#K2 = 0.001
#K4 = 0.00001
# Parâmetros do modelo (artigo de Mattos et al., 2021):
beta = 0.5
K2 = 0.1
K4 = 0.00001
# Constantes auxiliares (simplificação de notação):
kappa1 = beta*K2*dt/(dx**2)
kappa2 = - ((beta*(1-beta))*K4*dt)/(dx**2)
a1 = -kappa2
a2 = ((4*kappa2) - kappa1)
a3 = (1+(2*kappa1) - (6*kappa2))
# Vetor com as posições (para plotagem):
X = np.linspace(0.0, L, N)

totalsteps = int(tf3/dt)
pos025 = int(0.2*N)
pos050 = int(0.5*N)
pos075 = int(0.9*N)

tempos= np.linspace(0.0, tf3, totalsteps)
T1= np.linspace(0.0, tf3, totalsteps)
T2= np.linspace(0.0, tf3, totalsteps)
T3= np.linspace(0.0, tf3, totalsteps)

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

def plot_side_by_side(x, y1, y2):
    plt.subplot(1, 2, 1)
    plt.plot(x, y1, 'r')
    plt.title('Condição Inicial')
    plt.grid(True, 'major', 'both')
    plt.xlabel('x (m)')
    plt.ylabel('Concentração')

    plt.subplot(1, 2, 2)
    plt.plot(x, y2, 'b')
    plt.grid(True, 'major', 'both')
    plt.title('Solução analítica em t = 0')
    plt.xlabel('x (m)')

    plt.savefig('Comparados.png')
    plt.show()

def plot_evolve(x, y1a, y2a, y3a):
    plt.plot(x, y1a, 'blue', label='x = 0.2L', linewidth=2)
    plt.plot(x, y2a, 'red', label='x = 0.5L', linewidth=2)
    plt.plot(x, y3a, 'green', label='x = 0.9L', linewidth=2)
    plt.xlabel("Tempo (s)", fontsize = 11)
    plt.ylabel("Concentração", fontsize = 11)
    plt.grid(True, 'major', 'both')
    plt.legend(loc='best', fontsize=8)
    plt.savefig('Comparacao_pontos.png')
    plt.show() 

def plot_triple(x, y1a,y2a, y3a):
    plt.plot(x, y1a, 'blue', label='t = 0.5', linewidth=2)
    plt.plot(x, y2a, 'red', label='t = 4.0', linewidth=2)
    plt.plot(x, y3a, 'green', label='t = 10.0', linewidth=2)
    plt.xlabel("Posição (m)", fontsize = 11)
    plt.ylabel("Concentração", fontsize = 11)
    plt.grid(True, 'major', 'both')
    plt.legend(loc='best', fontsize=8)
    plt.savefig('Comparacao_tempos.png')
    plt.show() 
 
def plot_threepairs(x, y1a, y1b, y2a, y2b, y3a, y3b):
    plt.plot(x, y1a, 'blue', label='t = 10 (Analítico)', linewidth=2)
    plt.plot(x, y1b, 'red', label='t = 10 (Numérico)', linestyle='dashed', linewidth=2)
    plt.plot(x, y2a, 'purple', label='t = 50 (Analítico)', linewidth=2)
    plt.plot(x, y2b, 'lime', label='t = 50 (Numérico)', linestyle='dashed', linewidth=2)
    plt.plot(x, y3a, 'black', label='t = 100 (Analítico)', linewidth=2)
    plt.plot(x, y3b, 'orange', label='t = 100 (Numérico)', linestyle='dashed', linewidth=2)
    
    plt.xlabel("Posição (m)", fontsize = 11)
    plt.ylabel("Concentração", fontsize = 11)
    plt.grid(True, 'major', 'both')
    plt.legend(loc='best', fontsize=8)
    plt.savefig('Comparacao_Tripla.png')
    plt.show()

def initial_condition(x):
    return math.sin(math.pi*x/L)

def analytical_solution(x, t):
    return math.exp(-beta*K2*(math.pi**2)*t)*math.exp(-beta*(1-beta)*K4*(math.pi**4)*t)*math.sin(math.pi*x/L)

def solve_explicitly(tf):
    global T1; global T2; global T3;
    nsteps = int(tf/dt)
    print("Avaliando solução implicita para t = ", tf)
    print("Número de passos de tempo para calcular: ", nsteps)
 
    A =  np.zeros((N, N))         # Matriz dos coeficientes
    B = np.zeros(N)               # Matriz dos termos independentes

    # Preenchemos A:
    A[0, 0] = 1.0
    A[1, 0] = 2.0; A[1, 1] = -5.0; A[1, 2] = 4.0; A[1, 3] = -1.0;
    for i in range(2, N-2):
        A[i, i-2] = a1
        A[i, i-1] = a2
        A[i, i]   = a3
        A[i, i+1] = a2
        A[i, i+2] = a1
    A[N-2, N-4] = -1.0;  A[N-2, N-3] = 4.0; A[N-2, N-2] = -5.0; A[N-2, N-1] = 2.0;
    A[N-1, N-1] = 1.0

    # O vetor de termos independentes é, inicialmente, os P no instante 0, então preenchemos B com a condição inicial
    # Preenchemos da terceira (índice 2) até a antepenúltima (índice N-3) célula:
    for i in range (2, N-2):
        B[i] = initial_condition(i*dx)

    # Iteração no tempo:
    for n in range(0, nsteps):
        # Corrigimos os Bs:
        B[0] = 0.0; B[1] = 0.0; B[N-2] = 0.0; B[N-1] = 0.0
        # Obtemos o novo B como a solução do sistema linear:
        B = np.linalg.solve(A, B)
        T1[n] = B[pos025]
        T2[n] = B[pos050]
        T3[n] = B[pos075]
    return B

InitialCond = np.zeros(N)         # Condição inicial
ZeroSolution = np.zeros(N)        # Solução analítica em t = 0

for i in range (0, N):
    InitialCond[i] = initial_condition(i*dx)
    ZeroSolution[i] = analytical_solution(i*dx, 0.0)


Sol_at_10 = np.zeros(N)
Sol_at_50 = np.zeros(N)
Sol_at_100 = np.zeros(N)

# Solução analítica nos três tempos analisados
for i in range (0, N):
    Sol_at_10[i] = analytical_solution(i*dx, tf1)
    Sol_at_50[i] = analytical_solution(i*dx, tf2)
    Sol_at_100[i] = analytical_solution(i*dx, tf3)

# Solução numérica no tempo requerido:

Numeric_at_100 = solve_explicitly(tf3)
plot_evolve(tempos, T1, T2, T3)