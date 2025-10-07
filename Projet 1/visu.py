
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
from matplotlib import rc

# a, b = 1, 1

# f = open("U_sol.txt", "r")
# lines = f.readlines()

# Nx = int(np.sqrt(len(lines)) -2)    #on travail sur un carre
# Ny = Nx

# x = np.linspace(0, a, Nx+2)
# y = np.linspace(0, b, Ny+2)
# X, Y = np.meshgrid(x, y)

# U=np.zeros((Nx+2,Ny+2))

# for i in range(Nx+2):
#     for j in range(Ny+2):
#         U[i,j] = (float(lines[i*(Ny+2)+j]))


# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# ax.plot_surface(X, Y, U.T, cmap=cm.viridis)

# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("u(x,y)")


# xlim = (0,1)
# ylim = (0,1)
# zlim = (0,1)
# ax.set_xlim(xlim)
# ax.set_ylim(ylim)
# ax.set_zlim(zlim)
# ax.set(xticklabels=[],
#        yticklabels=[],
#        zticklabels=[])

# plt.show()

# ===============================================================================
# ===============================================================================

# f = open("norme_inf.txt", "r")
# lines = f.readlines()

# U=[]

# for line in lines: 
#     U.append(float(line))

# K = np.linspace(4,4-1+len(lines),len(lines))

# I = 2**K; I=1/I
# print(I)
# X2 = I**2
# X1 = I
# # X1 = 1/I; X2 = 1/X2
# # X1 /=U[0]; X2/=U[0]
# plt.loglog(I,U,'-o' ,label='err')
# plt.loglog(I,X1,'--',label='h')
# plt.loglog(I,X2,'--',label='h^2')
# plt.grid()
# plt.legend()
# plt.xlabel("dx")
# plt.ylabel("norme inf")

# plt.xlim(I[-1]/1.1,I[0]*1.1)
# plt.show()


# ===============================================================================
# ===============================================================================


f = open("iter_conv.txt", "r")
lines = f.readlines()

U=[]

for line in lines: 
    U.append(float(line))

K = np.linspace(4,4-1+len(lines),len(lines))
I = 2**K; I=1/I
X2 = I**2
X1 = I
X1 = 1/I; X2 = 1/X2

print(U,I)
plt.grid()
plt.loglog(I,U,'-o' ,label='iterations')
plt.loglog(I,X1,'--',label='h^-1')
plt.loglog(I,X2,'--',label='h^-2')
plt.legend()
plt.xlabel("dx")
plt.ylabel("norme inf")

plt.xlim(I[0]*1.1,I[-1]/1.1)
plt.show()

