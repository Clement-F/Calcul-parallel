import numpy as np
import matplotlib.pyplot as plt

# Ny = 64
# Nx = 32, 64, 96, 128, 160, 192, 224, 256
# P = 1, 2, 3, 4, 5, 6, 7, 8
ensta = True
cholesky = False

if ensta:
    P_ensta = 8
    P_list_ensta = np.arange(1,P_ensta+1)
    T_ensta = np.array([2.20469, 7.57985, 17.4631, 32.6088, 99.4945, 165.547, 253.299, 367.151])
    Efficiency_ensta = np.zeros(T_ensta.size)
    for i, t in enumerate(T_ensta):
        Efficiency_ensta[i] = T_ensta[0]/(P_list_ensta[i]*T_ensta[i])
    # Nx = 64, 128, 192, 256, 320, 384, 448, 512
    fig, ax = plt.subplots(figsize=(6,4), dpi=200)
    ax.plot(P_list_ensta, Efficiency_ensta, '-o', label='Mesuré', color='red')
    ax.plot(P_list_ensta, np.ones(P_ensta),'--' , label='Idéal', color='black')
    ax.legend()
    ax.grid()
    ax.set_xlabel('Nombre de Processeurs')
    ax.set_ylabel('Efficacité')
    plt.tight_layout()
    plt.show()

#T = np.araray([2.39765, 5.9221, 12.1449, 10.069, 66.1034, 39.19, 103.302, 77.6271])

# Cholesky
if cholesky:
    P_cholesky = 64
    P_list_cholesky = np.array([1, 2, 4, 8, 16, 24, 32, 48, 64])
    # Nx = 32, 64, 128, 256, 512, 768, 1024, 1536, 2048
    T_cholesky = np.array([2.08717, 3.41869, 9.37824, 43.7693, 154.773, 332.39, 369.724, 431.017, 422.52])
    Efficiency_cholesky = np.zeros(T_cholesky.size)
    for i, t in enumerate(T_cholesky):
        Efficiency_cholesky[i] = T_cholesky[0]/(P_list_cholesky[i]*T_cholesky[i])
    fig, ax = plt.subplots(figsize=(6,4), dpi=250)
    ax.plot(P_list_cholesky, Efficiency_cholesky, '-o', label='Mesuré', color='red')
    ax.plot(P_list_cholesky, np.ones(P_list_cholesky.size),'--' , label='Idéal', color='black')
    ax.legend()
    ax.grid()
    ax.set_xlabel('Nombre de Processeurs')
    ax.set_ylabel('Efficacité')
    plt.tight_layout()
    plt.show()
