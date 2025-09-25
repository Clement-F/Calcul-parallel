
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

Nx, Ny = 200, 200
a, b = 1, 1

x = np.linspace(0, a, Nx+2)
y = np.linspace(0, b, Ny+2)
X, Y = np.meshgrid(x, y)

u = np.loadtxt("u_sol.txt").reshape((Nx+2, Ny+2))


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_surface(X, Y, u.T, cmap=cm.viridis)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("u(x,y)")

ax.set(xticklabels=[],
       yticklabels=[],
       zticklabels=[])

plt.show()
