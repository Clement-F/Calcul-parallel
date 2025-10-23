
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
from matplotlib import rc

a, b = 1, 1

f = open("U_sol.txt", "r")
lines = f.readlines()

Nx = int(np.sqrt(len(lines)) -2)    #on travail sur un carre
Ny = Nx

x = np.linspace(0, a, Nx+2)
y = np.linspace(0, b, Ny+2)
X, Y = np.meshgrid(x, y)

U=np.zeros((Nx+2,Ny+2))

for i in range(Nx+2):
    for j in range(Ny+2):
        U[i,j] = (float(lines[i*(Ny+2)+j]))

m = np.min(U)
M = np.max(U)


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_surface(X, Y, U.T, cmap=cm.viridis)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("u(x,y)")


xlim = (0,1)
ylim = (0,1)
zlim = (m,M)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_zlim(zlim)
ax.set(xticklabels=[],
       yticklabels=[],
       zticklabels=[])

plt.show()
