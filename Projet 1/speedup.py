import numpy as np
import matplotlib.pyplot as plt


P = 8
P_list = np.arange(1,P+1)

t_seq = 37.2559
T_para = np.array([37.4174, 18.5754, 12.937, 10.053, 16.4673, 16.5504, 16.0822, 16.7312])
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