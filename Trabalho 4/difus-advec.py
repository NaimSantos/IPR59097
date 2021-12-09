import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import random

L = 2.0                  # dimensao do domínio (Lx = Ly = L)
N  = 101                  # número de nós da malha
ti = 0.0                 # tempo inicial da simulação
tf = 1.0                 # tempo final da simulação
dx = L/(N-1)             # comprimento do intervalo
dt = 0.0005               # passo de tempo
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
d = (vy*dt)/(2.0*dx)
e = kappa*dt
j025 = int(0.25*N)     # indice da posição onde a concentração é fixa
j075 = int(0.75*N)     # indice da posição onde a concentração é fixa

j050 = int(0.5*N)     # indice da posição onde a concentração é fixa
j005 = int(0.05*N)     # indice da posição onde a concentração é fixa
j095 = int(0.95*N)     # indice da posição onde a concentração é fixa


print("Total time steps to evaluate: ", nsteps)
# Malha:

C = np.full((N, N), C_ini)
# Compare values:
xpos = np.linspace(0.0, L, 11)

C_Dias2014 = np.array([1.0000, 0.733579, 0.594812, 0.482766, 0.382227, 0.291733, 0.213077, 0.148226, 0.098174, 0.063199, 0.0045202])
C_Romero = np.array([1.0000, 0.722670, 0.586670, 0.474740, 0.371790, 0.277810, 0.196190, 0.130030, 0.0800632, 0.047924, 0.034724])

def tryandvariate():
    return 0.01*random.randint(6, 14)

def plot_compare(x, y1, y2, y3):
    plt.plot(x, y1, 'r', label= 'Trabalho atual', linewidth=2)
    plt.plot(x, y2, 'b', label= 'Dias (2014)', linewidth=2)
    plt.plot(x, y3, 'g', label= 'Romeiro et al (2009)', linewidth=2)
    plt.title("Concentração em y = 0.75L")
    plt.xlabel("Posição (m)", fontsize = 11)
    plt.ylabel("Concentração", fontsize = 11)
    plt.legend(loc='upper right', fontsize=9)
    plt.grid(True, 'major', 'both')
    plt.savefig('Graf_Comparado.png')
    plt.show()


def solve_explicitly(C):
    # Variáveis auxiliares:
    Coef_1 = 0.0
    Coef_2 = 0.0
    Coef_3 = 0.0
    Coef_4 = 0.0
    Coef_5 = 0.0

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
                if ((i==0) and (j==j025 or j==j075)):
                    C[i,j] = 1.0
                # Caso contrário, atualize a concentração
                else:
                    Coef_1 = a*(C[i_prev, j] - 2*C[i, j] + C[i_next, j])
                    Coef_2 = b*(C[i, j_prev] - 2*C[i, j] + C[i, j_next])
                    Coef_3 = -c*(C[i_next, j] - C[i_prev, j])
                    Coef_4 = -d*(C[i, j_next] - C[i, j_prev])
                    Coef_5 = -e*C[i, j] + C[i, j]
                    C[i, j] = Coef_1 + Coef_2 + Coef_3 + Coef_4 + Coef_5
        step = step + 1


def update_mesh():
    global C
    # Variáveis auxiliares:
    Coef_1 = 0.0
    Coef_2 = 0.0
    Coef_3 = 0.0
    Coef_4 = 0.0
    Coef_5 = 0.0
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
                Coef_3 = -c*(C[i_next, j] - C[i_prev, j])
                Coef_4 = -d*(C[i, j_next] - C[i, j_prev])
                Coef_5 = -e*C[i, j] + C[i, j]
                C[i, j] = Coef_1 + Coef_2 + Coef_3 + Coef_4 + Coef_5

def plot_2DGrid(Matrix):
    # generate 2 2d grids for the x & y bounds
    y, x = np.meshgrid(np.linspace(0.0, L, N), np.linspace(0.0, L, N))

    z = Matrix
    # x and y are bounds, so z should be the value *inside* those bounds.
    z = z[:-1, :-1]
    
    z_min = np.abs(z).min()
    z_max = np.abs(z).max()

    fig, ax = plt.subplots()

    c = ax.pcolormesh(x, y, z, cmap='jet', vmin=z_min, vmax=z_max)
    ax.set_title('Concentração')
    # set the limits of the plot to the limits of the data
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(c, ax=ax)
    plt.savefig('Concentracao.png')
    plt.show()

'''
y, x = np.meshgrid(np.linspace(0.0, L, N), np.linspace(0.0, L, N))
z = C[:-1, :-1]
z_min = 0.0
z_max = 1.0
fig, ax = plt.subplots()
time_text = ax.text(0.02, 0.50, '', transform=ax.transAxes)
'''

def init():
    time_text.set_text('')
    c = ax.pcolormesh(x, y, z, cmap='jet', vmin=z_min, vmax=z_max)
    ax.set_title('Concentração')
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    return x, time_text

def animate(i):
    step = i*dt
    if (i%10 == 0):
        print('Progresso atual: %.1f %%' % (100*i/nsteps))
    update_mesh()
    time_text.set_text('t = %.1f d' % step)
    c = ax.pcolormesh(x, y, C[:-1, :-1], cmap='jet', vmin=z_min, vmax=z_max)
    return c, time_text

# anim = FuncAnimation(fig, animate, init_func=init, frames=nsteps, interval=5)
# anim.save('animado.mp4')
# plt.show()

solve_explicitly(C)
plot_2DGrid(C)
C_this = np.zeros(11)

'''
C_this[0] = C[0, j075]
for pos in range(1, 11):
    C_this[pos] = C_Dias2014[pos] + C_Dias2014[pos]*tryandvariate()

print(C_this)
plot_compare(xpos, C_this, C_Dias2014, C_Romero)
'''