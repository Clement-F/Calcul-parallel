
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

Ny=100; Nx=100
a=1;    b=1

X = np.linspace(0, a, Nx+1)
Y = np.linspace(0, b, Ny+1)
X, Y = np.meshgrid(X, Y)

f = open("U_sol.txt", "r")
lines = f.readlines()

U=np.zeros((Nx+1,Ny+1))

for i in range(Nx+1):
    for j in range(Ny+1):
        U[i,j] = (float(lines[i*Ny+j]))


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_surface(X, Y, U, vmin=U.min() * 2)

ax.set(xticklabels=[],
       yticklabels=[],
       zticklabels=[])

plt.show()