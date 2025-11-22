import numpy as np
import matplotlib.pyplot as plt

convergence = True
iterations = True

if convergence:
    f = open("tp1-exo2-convergence_0025.txt", "r"); lines_1 = f.readlines()
    f = open("tp1-exo2-convergence_001.txt", "r"); lines_2 = f.readlines()
    f = open("tp1-exo2-convergence_0005.txt", "r"); lines_3 = f.readlines()
    f = open("tp1-exo2-convergence_00025.txt", "r"); lines_4 = f.readlines()

    data1 = []; data2 = []; data3 = []; data4 = []
    for l in lines_1: data1.append(float(l))
    for l in lines_2: data2.append(float(l))
    for l in lines_3: data3.append(float(l))
    for l in lines_4: data4.append(float(l))

    for data in [data1, data2, data3, data4]:
        for i in range(len(data)):
            if data[i] != 0:
                dernier_nonnul = data[i]
            else:
                data[i] = dernier_nonnul


    k = 250
    k_list = np.arange(1,k)

    fig, ax = plt.subplots(figsize=(6,4), dpi=200)
    ax.plot(k_list, data1[:k], label='$h=0.01$')
    ax.plot(k_list, data2[:k], label='$h=0.025$')
    ax.plot(k_list, data3[:k], label='$h=0.0025$')
    ax.plot(k_list, data4[:k], label='$h=0.0005$')
    ax.grid()
    ax.legend()
    ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.set_xlabel("Nombre d'itérations")
    ax.set_ylabel('$\\frac{|| u_h^{(k)} - u_h ||_{H^1(\\Omega)}}{||u_h||_{H^1(\\Omega)}}$')
    ax.set_title("Erreur d'approximations par la méthode du gradient \n conjugué en fonction du nombre d'itérations")
    plt.tight_layout()
    plt.show()


if iterations:
    data = np.loadtxt("tp1-exo2-iterations.txt")
    h = data[:, 0]
    iter = data[:, 1]
    fig, ax = plt.subplots(figsize=(6,4), dpi=250)
    ax.plot(h, iter)
    ax.grid()
    ax.set_xlabel('$h$')
    ax.set_ylabel("Nombre d'itérations de CG nécessaires pour \n atteindre $\\frac{|| u_h^{(k)} - u_h ||_{H^1(\\Omega)}}{||u_h||_{H^1(\\Omega)}} < 10^{-6}$")
    plt.tight_layout()
    plt.show()