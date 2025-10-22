import numpy as np
import matplotlib.pyplot as plt

# Ny = 64
# Nx = 32, 64, 96, 128, 160, 192, 224, 256
# P = 1, 2, 3, 4, 5, 6, 7, 8
ensta = True
cholesky = True

if ensta:
    P_ensta = 8
    P_list_ensta = np.arange(1,P_ensta+1)
    T_ensta = np.array([1.9861, 2.08667, 2.23234, 2.36097, 4.53718, 5.04197, 5.82658, 7.03631])
    Efficiency_ensta = np.zeros(T_ensta.size)
    for i, t in enumerate(T_ensta):
        Efficiency_ensta[i] = T_ensta[0]/(T_ensta[i])
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

# Cholesky
if cholesky:
    P_cholesky = 64
    P_list_cholesky = np.array([1, 2, 4, 8, 16, 24, 32, 48, 64])
    # Nx = 32, 64, 128, 256, 512, 768, 1024, 1536, 2048
    T_cholesky = np.array([2.47246, 2.53846, 2.84465, 2.97019, 3.0464, 3.10246, 3.13227, 3.22837, 3.29659])
    Efficiency_cholesky = np.zeros(T_cholesky.size)
    for i, t in enumerate(T_cholesky):
        Efficiency_cholesky[i] = T_cholesky[0]/(T_cholesky[i])
    fig, ax = plt.subplots(figsize=(6,4), dpi=250)
    ax.plot(P_list_cholesky, Efficiency_cholesky, '-o', label='Mesuré', color='red')
    ax.plot(P_list_cholesky, np.ones(P_list_cholesky.size),'--' , label='Idéal', color='black')
    ax.legend()
    ax.grid()
    ax.set_xlabel('Nombre de Processeurs')
    ax.set_ylabel('Efficacité')
    plt.tight_layout()
    plt.show()
