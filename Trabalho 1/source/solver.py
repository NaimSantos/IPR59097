import math
import numpy as np
import matplotlib.pyplot as plt

kappa = 59.0                        # coeficiente de condutividade térmica
h = 380                             # coeficiente de troca de calor por convecção
D = 4e-03                           # diâmetro da seção transversal da aleta
T_amb = 20.0                        # temperatura do ambiente
T_0 = 320.0                         # temperatura da aleta em x=0
T_n = 25.0                          # temperatura da aleta em x=L
L = 0.25                            # comprimento da aleta
N = 30                              # número de nós em x
dx = L/(N-1)                        # comprimento do intervalo em x
area = math.pi*D*D/4                # área da seção transversal da aleta
peri = math.pi*D                    # perímetro da aleta
m2 = (h*peri)/(kappa*area)
gamma = dx*dx*m2
m = math.sqrt(m2)

# Matrizes para armazenamento dos resultados
A = np.zeros((N, N))                # Array para os coeficientes, com NxN elementos
B = np.zeros(N)                     # Vetor dos termos independentes
X = np.zeros(N)                     # Vetor para armazenar as temperaturas

# Corrige a matriz de coeficientes:
A[0][0] = 1.0
i = 1
while i < (N-1) :
    A[i][i-1] = 1.0
    A[i][i] = - (2 + gamma)
    A[i][i+1] = 1.0
    i = i + 1
A[N-1][N-1] = 1.0
#print(A)


# Corrige os termos independentes:
B[0] = T_0
j = 1
while j < N-1 :
    B[j] = -gamma * T_amb
    j = j + 1
B[N-1] = T_n
#print(B)

# Armazena em X a solução do sistema
X = np.linalg.solve(A, B)
print(X)
Pos = np.linspace(0.0, L, N)        # Vetor das posições linearmente espaçado


plt.plot(Pos, X, 'r')
plt.title("Perfil de temperatura da placa unidimensional (t=500 s)")
plt.xlabel("Comprimento (m)", fontsize = 11)
plt.ylabel("Temperatura (° C)", fontsize = 11)
plt.grid(True, 'major', 'both')
plt.savefig('temperatura.png')
plt.show()

