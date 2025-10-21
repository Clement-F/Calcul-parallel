import numpy as np
import matplotlib.pyplot as plt


ensta = False

if ensta:
    # Nx = Ny = 128
    P = 8
    P_list = np.arange(1,P+1)

    t_seq = 37.2559
    T_para = np.array([37.4174, 18.5754, 12.937, 10.053, 16.4673, 16.5504, 16.0822, 16.7312])
    
    Speedup = np.zeros(T_para.size)
    for i, t in enumerate(T_para):
        Speedup[i] = t_seq/t
else:
    # Nx = Ny = 256
    P = 64
    P_list = np.array([1, 2, 3, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64])
    t_seq = 1582.67
    T_para = np.array([1582.67, 1071.86, 735.855, 576.587, 299.574, 162.446, 128.234, 103.034,  85.9106, 79.0844, 74.045, 80.0096, 79.1271, 75.591, 67.6915, 70.4882, 64.1522, 68.001, 61.5305])
    Speedup = np.zeros(T_para.size)
    for i, t in enumerate(T_para):
        Speedup[i] = t_seq/t



fig, ax = plt.subplots(figsize=(6,4), dpi=150)
ax.plot(P_list, P_list, '--', label='Idéal', color='black')
ax.plot(P_list, Speedup, label='Mesuré', color='red')
ax.set_xlabel('Nombre de Processeurs')
ax.set_ylabel('Speedup')
ax.legend()
ax.set_xlim(1, P)
ax.set_ylim(1, P)
ax.grid()
plt.tight_layout()
plt.show()
