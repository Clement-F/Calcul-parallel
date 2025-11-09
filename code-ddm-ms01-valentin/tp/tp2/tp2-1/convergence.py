import numpy as np
import matplotlib.pyplot as plt

convergence = True
iterations = True

if convergence:
    f = open("convergence_001.txt", "r"); lines_1PCG = f.readlines()
    f = open("convergence_0025.txt", "r"); lines_2PCG = f.readlines()
    f = open("convergence_0005.txt", "r"); lines_3PCG = f.readlines()
    f = open("convergence_00025.txt", "r"); lines_4PCG = f.readlines()

    data1 = []; data2 = []; data3 = []; data4 = []
    for l in lines_1PCG: data1.append(float(l))
    for l in lines_2PCG: data2.append(float(l))
    for l in lines_3PCG: data3.append(float(l))
    for l in lines_4PCG: data4.append(float(l))

    for data in [data1, data2, data3, data4]:
        for i in range(len(data)):
            if data[i] != 0:
                dernier_nonnul = data[i]
            else:
                data[i] = dernier_nonnul

    k = 100
    k_list = np.arange(1,k+1)

    fig, ax = plt.subplots(figsize=(6,4), dpi=200)
    ax.plot(k_list, data1[:k], label='$h=0.01$')
    ax.plot(k_list, data2[:k], label='$h=0.025$')
    ax.plot(k_list, data3[:k], label='$h=0.0025$')
    ax.plot(k_list, data4[:k], label='$h=0.0005$')
    ax.grid()
    ax.legend()
    ax.set_yscale('log')
    ax.set_xlabel("Nombre d'itérations")
    ax.set_ylabel('$\\frac{|| u_h^{(k)} - u_h ||_{H^1(\\Omega)}}{||u_h||_{H^1(\\Omega)}}$')
    ax.set_title("Erreur d'approximations par la méthode du gradient \n conjugué préconditionné en fonction du nombre d'itérations")
    plt.tight_layout()
    plt.show()

if iterations:
    iter_CG = np.loadtxt("../../tp1/tp1-3/iterations.txt")
    iter_PCG = np.loadtxt("iterations_PCG.txt")
    h = iter_PCG[:, 0]
    iterations_CG = iter_CG[:, 1]
    iterations_PCG = iter_PCG[:, 1]
    fig, ax = plt.subplots(figsize=(6,4), dpi=250)
    ax.plot(h, iterations_PCG, '-o', label='PCG')
    ax.plot(h, iterations_CG, '-o', label='CG')
    ax.grid()
    ax.set_xlabel('$h$')
    ax.set_ylabel("Nombre d'itérations nécessaires pour \n atteindre $\\frac{|| u_h^{(k)} - u_h ||_{H^1(\\Omega)}}{||u_h||_{H^1(\\Omega)}} < 10^{-6}$")
    ax.legend()
    plt.tight_layout()
    plt.show()