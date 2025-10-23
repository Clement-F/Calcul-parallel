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

    #mesures
    T_ensta_J = np.array([1.9861, 2.08667, 2.23234, 2.36097, 4.53718, 5.04197, 5.82658, 7.03631])
    T_ensta_GS = np.array([6.48086, 13.5267, 21.8505, 30.1659, 72.3862, 90.7497, 108.626, 124.006])
    Efficiency_ensta_J = np.zeros(P_list_ensta.size)
    Efficiency_ensta_GS = np.zeros(P_list_ensta.size)

    for i, t in enumerate(T_ensta_J):
        Efficiency_ensta_J[i] = T_ensta_J[0]/(T_ensta_J[i])
    
    
    for i, t in enumerate(T_ensta_J):
        Efficiency_ensta_GS[i] = T_ensta_GS[0]/(T_ensta_GS[i])

    # Nx = 64, 128, 192, 256, 320, 384, 448, 512
    fig, ax = plt.subplots(figsize=(6,4), dpi=200)
    ax.plot(P_list_ensta, Efficiency_ensta_J,  '-o', label='Jacobi', color='green')
    ax.plot(P_list_ensta, Efficiency_ensta_GS, '-o', label='Gauss-Seidel', color='red')
    ax.plot(P_list_ensta, np.ones(P_ensta),'--' , label='Idéal', color='black')
    ax.legend()
    ax.grid()
    ax.set_xlabel('Nombre de Processeurs')
    ax.set_ylabel('Efficacité')
    ax.set_ylim(0,1.1); ax.set_xlim(1,P_ensta)
    plt.tight_layout()
    plt.savefig("Comp_Efficency_Ensta.png")

# Cholesky
if cholesky:
    P_cholesky = 64
    P_list_cholesky = np.array([1, 2, 4, 8, 16, 24, 32, 48, 64])
    # Nx = 32, 64, 128, 256, 512, 768, 1024, 1536, 2048
    # mesures
    T_cholesky_J = np.array([2.47246, 2.53846, 2.84465, 2.97019, 3.0464, 3.10246, 3.13227, 3.22837, 3.29659])
    T_cholesky_GS = np.array([ 6.12499,  13.6001, 31.1099, 70.2722, 121.313, 190.981, 250.104, 434.25, 586.501])  

    Efficiency_cholesky_J = np.zeros(P_list_cholesky.size)
    Efficiency_cholesky_GS = np.zeros(P_list_cholesky.size)

    for i, t in enumerate(T_cholesky_J):
        Efficiency_cholesky_J[i] = T_cholesky_J[0]/(T_cholesky_J[i])
        
    for i, t in enumerate(T_cholesky_GS):
        Efficiency_cholesky_GS[i] = T_cholesky_GS[0]/(T_cholesky_GS[i])

    fig, ax = plt.subplots(figsize=(6,4), dpi=250)
    ax.plot(P_list_cholesky, Efficiency_cholesky_J, '-o', label='Jacobi', color='green')
    ax.plot(P_list_cholesky, Efficiency_cholesky_GS, '-o', label='Gauss-Seidel', color='red')
    ax.plot(P_list_cholesky, np.ones(P_list_cholesky.size),'--' , label='Idéal', color='black')
    ax.legend()
    ax.grid()
    ax.set_xlabel('Nombre de Processeurs')
    ax.set_ylabel('Efficacité')
    ax.set_ylim(0,1.1); ax.set_xlim(1,P_cholesky)
    plt.tight_layout()
    plt.savefig("Comp_Efficency_Chol.png")
