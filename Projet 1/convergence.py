import numpy as np
import matplotlib.pyplot as plt

a, b = 1, 1
noPara = True
GaussSeidel = True
Jacobi = False

N_list = np.array([8, 16, 32, 64, 128])
h_list = np.zeros(N_list.size)
for i, N in enumerate(N_list):
    h_list[i] = a/(N_list[i]+1)

conver_iter = np.array([950, 3596, 13968, 55035])
if Jacobi:
    if noPara:
        # Erreurs sans parallélisation
        err_list = np.array([0.00990807, 0.00282655, 0.000753932, 0.000199798, 4.94689e-05])
    else:
        # Erreurs avec parallélisation, P = 1
        err_list = np.array([0.00990807, 0.00282655, 0.000753932, 0.000199798, 4.94689e-05])

if GaussSeidel:
    # Erreurs Gauss-Seidel
    err_list2 = np.array([0.00990797, 0.00282644, 0.000753813, 0.000194506, 4.93482e-05])


fig, ax = plt.subplots(figsize=(6,4), dpi=250)

ax.plot(h_list, h_list**2, '--', label="$h^2$")
ax.plot(h_list, h_list, '--', label="$h$")
if Jacobi:
    ax.plot(h_list, err_list, '-o', label='Erreur')

if GaussSeidel:
    ax.plot(h_list, err_list2, '-o', label='Erreur')

ax.set_xlabel('$h = max(\Delta x, \Delta y)$')
ax.set_ylabel('Erreur $L^\infty$')
ax.loglog()
ax.grid()
ax.legend()
plt.tight_layout()
plt.show()