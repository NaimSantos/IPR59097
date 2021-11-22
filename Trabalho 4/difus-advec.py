import numpy as np
import matplotlib.pyplot as plt
import math


L = 2.0                  # dimensao do domínio (Lx = Ly = L)
N  = 41                  # número de nós da malha
ti = 0.0                 # tempo inicial da simulação
tf = 1.0                 # tempo final da simulação
dx = L/(N-1)             # comprimento do intervalo
dt = 0.0005              # passo de tempo
nsteps = int((tf-ti)/dt) # número de passos de tempo
C_ini = 0.0              # concentração inicial
Dxx = 0.2
Dyy = 0.2
kappa = 0.001
vx = 1.0
vy = 0.0
a = (dt*Dxx)/(dx*dx)
b = (dt*Dyy)/(dx*dx)
c = (vx*dt)/(2.0*dx)
d = kappa*dt
j025 = int(0.25*N)     # indice da posição onde a concentração é fixa
j075 = int(0.75*N)     # indice da posição onde a concentração é fixa

# Malha:

C = np.full((N, N), C_ini)

def show_parameters():
    print('Refimamento: %.2f (m)' %dx)
    print("Tempo total de simulação: ", tf)
    print("Passo de tempo: ", dt)
    print("Número de passos de tempo: ", nsteps)


def plot_2DGrid(Matrix):
    # generate 2 2d grids for the x & y bounds
    y, x = np.meshgrid(np.linspace(0.0, L, N), np.linspace(0.0, L, N))

    z = Matrix
    # x and y are bounds, so z should be the value *inside* those bounds.
    z = z[:-1, :-1]
    
    z_min = np.abs(z).min()
    z_max = np.abs(z).max()

    fig, ax = plt.subplots()

    c = ax.pcolormesh(x, y, z, cmap='bwr', vmin=z_min, vmax=z_max)
    ax.set_title('Concentração')
    # set the limits of the plot to the limits of the data
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(c, ax=ax)
    plt.savefig('Concentracao.png')
    plt.show()


def solve_explicitly(C):
    # Variáveis auxiliares:
    Coef_1 = 0.0
    Coef_2 = 0.0
    Coef_3 = 0.0
    Coef_4 = 0.0

    # Evolução temporal:
    step = 0
    while (step < nsteps):

        if (step%100 == 0):
            print("Passo de tempo atual: ", step)

        # Iteração na malha
        for i in range (0, N):
            for j in range (0, N):

                # Implementação da periodicidade
                if (i==0):
                   i_prev = 0 # periodicidade real: i_prev = N-1
                else:
                    i_prev = i-1
                if (i==N-1):
                    i_next = N-1 # periodicidade real: i_next = 0
                else:
                    i_next = i+1
                # periodicidade em y
                if (j==0):
                    j_prev = N-1
                else:
                    j_prev = j-1
                if (j==N-1):
                    j_next = 0
                else:
                    j_next = j+1

                # Neste ponto, a concentração é dada:
                if (i==0 and (j==j025 or j==j075)):
                    C[i,j] = 1.0
                # Caso contrário, atualize a concentração
                else:
                    Coef_1 = a*(C[i_prev, j] - 2*C[i, j] + C[i_next, j])
                    Coef_2 = b*(C[i, j_prev] - 2*C[i, j] + C[i, j_next])
                    Coef_3 = -c*(C[i_prev, j] - C[i_next, j])
                    Coef_4 = d*C[i, j] + C[i, j]
                    C[i, j] = Coef_1 + Coef_2 + Coef_3 + Coef_4
        step = step + 1



show_parameters()
solve_explicitly(C)
plot_2DGrid(C)
