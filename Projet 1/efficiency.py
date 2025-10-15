import numpy as np
import matplotlib.pyplot as plt


P = 8
P_list = np.arange(1,P+1)

T = np.array([2.39765, 5.9221, 12.1449, 10.069, 66.1034, 39.19, 103.302, 77.6271])
Efficiency = np.zeros(T.size)

for i, t in enumerate(T):
    Efficiency[i] = T[0]/T[i]

print(Efficiency.size)

fig, ax = plt.subplots(figsize=(6,4), dpi=250)
ax.plot(P_list, Efficiency, label='Mesuré', color='red')
ax.plot(P_list, np.ones(P),'--' , label='Idéal', color='black')
ax.legend()
ax.grid()
ax.set_xlabel('Nombre de Processeurs')
ax.set_ylabel('Efficacité')
plt.tight_layout()
plt.show()