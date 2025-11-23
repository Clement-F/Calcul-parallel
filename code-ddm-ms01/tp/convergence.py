
import numpy as np
import matplotlib.pyplot as plt

f = open("convergence_0025.txt", "r"); lines_1 = f.readlines()
f = open("convergence_001.txt",  "r"); lines_2 = f.readlines()
f = open("convergence_0005.txt", "r"); lines_3 = f.readlines()
f = open("convergence_00025.txt","r"); lines_4 = f.readlines()

data1=[];data2=[];data3=[];data4=[]
for l in lines_1: data1.append(float(l))
for l in lines_2: data2.append(float(l))
for l in lines_3: data3.append(float(l))
for l in lines_4: data4.append(float(l))

h_list = np.arange(1,len(lines_1)+1)

fig, ax = plt.subplots(figsize=(6,4), dpi=250)

# ax.plot(h_list, 1/h_list**4, '--', label="$k^{-4}$")
# ax.plot(h_list, 1/h_list**3, '--', label="$k^{-3}$")
# ax.plot(h_list, 1/h_list, '--', label="$k^{-1}$")


ax.plot(h_list, data1, '-o', label='h= 0.025')
ax.plot(h_list, data2, '-o', label='h= 0.01')
ax.plot(h_list, data3, '-o', label='h= 0.005')
ax.plot(h_list, data4, '-o', label='h= 0.0025')

ax.set_xlabel('$k-iterations$')
ax.set_ylabel('Erreur $L^\infty$')
ax.loglog()
ax.grid()
ax.legend()
plt.tight_layout()
plt.savefig("convergence.png")

