
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm

a, b = 1, 1

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


f = open("norme_inf.txt", "r")
lines = f.readlines()

U=[]

for line in lines: 
    U.append(float(line))

K = np.arange(len(lines))
I = 2**K

plt.plot(np.log(I),np.log(U))
plt.show()


# def plot_film(psi, duration=10, frames_per_second=30, a=1,b=1):
    
#     fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#     t_data = np.linspace(0, 1, np.size(psi, 1)) # 1 is arbitrary here
#     x_data = np.linspace(0, a, Nx+1, endpoint=False)
#     y_data = np.linspace(0, b, Ny+1, endpoint=False)
#     X, Y = np.meshgrid(x_data, y_data)
#     # set the min and maximum values of the plot, to scale the axis
#     m = min(0, np.min(np.real(psi)), np.min(np.imag(psi)))
#     M = np.max(np.abs(psi))
    
#     # set the axis once and for all
#     ax.set(xlim=[0,a],ylim=[0,b], zlim=[m,M], xlabel='x', ylabel='psi')

#     # dummy plots, to update during the animation
#     sol_plot =ax.plot_surface(X,Y, U.T, cmap=cm.viridis)
#     ax.set_xlabel("x")
#     ax.set_ylabel("y")
#     ax.set_zlabel("u(x,y)")

#     # define update function as an internal function (that can access the variables defined before)
#     # will be called with frame=0...(duration*frames_per_second)-1
#     def update(frame):
#         print(frame)

#         # get the data by linear interpolation
#         # t = frame / (duration * frames_per_second)
#         # psi_t = np.array([np.interp(t, t_data, psi[i,j, :]) for i in range(Nx+1) for j in range(Ny+1)])

#         k+=1
#         psi_t = psi[:,:,k]

#         # update the plots
#         sol_plot.set_ydata(psi_t)

#     ani = animation.FuncAnimation(fig=fig, func=update, frames=duration*frames_per_second, interval=1000/frames_per_second)
#     return ani


# U=np.zeros((Nx+1,Ny+1))
# U_film =np.zeros((Nx+1,Ny+1,10))

# for k in range(10):
    
#     f = open("U_sol.txt", "r")
#     lines = f.readlines()

#     for i in range(Nx+1):
#         for j in range(Ny+1):
#             U[i,j] = (float(lines[i*Ny+j]))
#     U_film[:,:,k] =U

# ani =plot_film(U_film)
# plt.show()
