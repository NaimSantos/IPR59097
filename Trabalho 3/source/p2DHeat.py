import numpy as np
import matplotlib.pyplot as plt
import math

# Variáveis do problema e do domínio:
L = 0.1                  # comprimento total da placa ( Lx = Ly)
N = 25                   # número de nós da malha (considerando N = M)
tf = 120                # tempo final da simulação
dx = L/(N-1)             # comprimento do intervalo (dx = dy)
dt = 0.2                 # passo de tempo
nsteps = int(tf/dt)      # número de passos de tempo
kappa = 0.6
rho = 600
cp = 1200
h = 15.0
g = 100000
T_ini = 20.0             # temperatura inicial da placa
T_f = 20.0               # temperatura do fluido na área vazada
T1 = 100.0               # temperatura nas arestas superior e direita
alpha = kappa/(rho*cp)
r = (kappa*dt)/(rho*cp*dx*dx)
lambda = (alpha*g*dt)/kappa
gamma = 2*h*dx/kappa


# Malha e outros arrays:
tim = np.linspace(0.0, tf, nsteps)   # Tempos
Temp = np.full((N, N), T_ini)           # Malha NxN de temperaturas, inicializadas em T_ini.
test_array = np.arange(100*100).reshape(100,100)


# Definições de funções utilizadas:
def heatmap2d(arr: np.ndarray):
    plt.imshow(arr, cmap='YlOrRd')
    plt.colorbar()
    plt.savefig('Heat_Map.png')
    plt.show()

def show_parameters():
    print('\nDomínio: %f (m)' %L)
    print('Refimamento: %.2f (m)' %dx)
    print("Critério de estabilidade: ", r)
    print('Tempo total de simulação: %0.1f (s)' %tf)
    print('Passo de tempo: %.1f (s)' %dt)
    print("Número de passos de tempo: ", nsteps)

def set_outer_boundaries(Temp):
    # Aresta esquerda:
    for i in range(0, N-1):
        j = 0
        Temp[i, j] = T1
        #Temp[i, j] = (1/3)*(4*Temp[i, j+1] - Temp[i, j+2])

    # Aresta inferior:
    for j in range(1, N-1):
        i = 0
        Temp[i, j] = T1
        #Temp[i, j] = (1/3)*(4*Temp[i+1, j] - Temp[i+2, j])

    # Aresta superior:
    for i in range (0, N):
        j = N-1
        Temp[i, j] = T1

    # Aresta direita:
    for j in range(0, N):
        i = N-1
        Temp[i, j] = T1
 



def solve_explicitly(Temp):
    p = int(0.3*N)
    p_w = int(0.7*N)

    # Iteração no tempo:
    t_current = 0.0
    while (t_current < tf):

        # impomos as condições nos contornos externos:
        set_outer_boundaries(Temp)

        # Iteração em x, de 1 a N-2:
        for i in range (1, N-1):

            # Iteração em y (de 1 a N-2):
            for j in range (1, N-1):

                # Temperatura na região vazada
                if ((i>p and i<p_w) and (j>p and j<p_w)):
                    Temp[i, j] = T_f

                # Temperatura nas regiões do contorno interno (área vazada)
                elif (i==p+1 and (j>p and j<p_w)):
                    Temp[i, j] = (1/(3 + gamma))*(gamma*T_ini + 4*Temp[i-1, j] - Temp[i-2, j] )

                elif (i==p_w-1 and (j>p and j<p_w)):
                    Temp[i, j] = (1/(gamma - 3))*(gamma*T_ini - 4*Temp[i+1, j] + Temp[i+2, j] )

                elif (j==p+1 and (i>=p and i<=p_w)):
                    Temp[i, j] = (1/(3 + gamma))*(gamma * T_ini + 4*Temp[i, j-1] - Temp[i, j-2] )

                elif(j==p_w-1 and (i>=p and i<=p_w)):
                    Temp[i, j] = (1/(gamma - 3))*(gamma*T_ini - 4*Temp[i, j+1] + Temp[i, j+2])

                # Temperatura nos demais nós
                else:
                    Temp[i, j] = Temp[i, j] + r*(Temp[i+1, j] + Temp[i-1, j] + Temp[i, j+1] + Temp[i, j-1] - 4*Temp[i, j]) + lambda

        # Próximo passo de tempo
        t_current = t_current + dt




# Chamada das rotinas para resolver o problema:
solve_explicitly(Temp)
#Temp = np.transpose(Temp)
heatmap2d(Temp)



