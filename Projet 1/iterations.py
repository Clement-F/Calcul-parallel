import numpy as np
import matplotlib.pyplot as plt

N_list = np.array([8, 16, 32, 64, 128])
h_list = np.zeros(N_list.size)
for i, N in enumerate(N_list):
    h_list[i] = 1/(N_list[i]+1)


iter_list_J = np.array([262, 950, 3596, 13968, 55035])
iter_list_GS = np.array([135, 482, 1817, 7051, 27774])

fig, ax = plt.subplots(figsize=(6,4),dpi=250)


ax.plot(N_list, iter_list_J, '-o', label='Jacobi')
ax.plot(N_list, iter_list_GS, '-o', label='Gauss-Seidel')
ax.plot(N_list, N_list**2, '--', label='$N_x^2$')
ax.plot(N_list, N_list, '--', label='$N_x$')

ax.set_xlabel('$N_x$')
ax.set_ylabel("Nombre d'itérations")
ax.legend()
ax.grid()
ax.loglog()
plt.tight_layout()
plt.show()

