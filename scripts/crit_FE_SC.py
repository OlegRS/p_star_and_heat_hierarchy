#####################################################
# This script plots the free energy at the critical #
# point alongside the critical isotherm.            #
#####################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

Ld = .3
fig, axs = plt.subplots(2,1, figsize=(10., 6.*2))

t2_MF = 1.8707482993197284
t3_MF = -0.7558578987150416
t1_MF = 1/2*np.log(Ld/(1-Ld)) - 2*t2_MF*Ld - 3*t3_MF*Ld*Ld


title = "$t_1=$" + str(round(t1_MF, 2)) + "$;\ t_2=$" + str(round(t2_MF,3)) + "$;\ t_3=$" + str(round(t3_MF,3))

L = np.arange(0,1, 1e-5)

F = -2*t1_MF*L - 2*t2_MF*L**2 - 2*t3_MF*L**3 + L*np.log(L/(1-L)) + np.log(1-L)
axs[0].plot(L, F, label="Free energy", linewidth=2, zorder=10, alpha=.8)
axs[0].plot(Ld, F[30000], marker='+', markersize=15, color='red', zorder=500, markeredgewidth=2, label="Critical point", linestyle='none')
axs[0].text(Ld+.005, F[30000]+.015, "(0.3, -0.1016)", size=16, color='r')

axs[0].axvline(0, color='black')
axs[0].axhline(0, color='black')

axs[0].set_xlim([0, 1])
axs[0].set_ylim([-.15, .5])
axs[0].axvline(0, color="black")
axs[0].axhline(0, color="black")
axs[0].legend(loc='upper left', prop={'size': 21})

# axs[0].title(title, fontsize=22)
axs[0].set_xlabel(r'$\langle L\rangle$', fontsize=26, labelpad=-8)
axs[0].set_ylabel(r'$\Phi(\langle L\rangle|{\bf t})$', fontsize=26, labelpad=-9)

plt.sca(axs[0])
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)

###########################################################################################
n_nodes = 1000

title = r"$n=10^3$" + "$;\ t_2=$" + str(round(t2_MF,3)) + "$;\ t_3=$" + str(round(t3_MF,3))

x_min, x_max =-2, .1
delta_t1 = 0.05
delta_L = 0.00001
t1_range = np.arange(x_min, x_max+delta_t1, delta_t1)
L_range = np.arange(0, 1+5*delta_L, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 +2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = axs[1].contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=1.5)
CS.collections[0].set_label('MF theor')

axs[1].set_xlabel(r'$t_1$', fontsize=26, labelpad=-1)
axs[1].set_ylabel(r'$\langle L\rangle$', fontsize=26, labelpad=-8)
axs[1].axvline(0, color="black")
axs[1].axhline(0, color="black")
axs[1].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
axs[1].set_ylim([-.01, 1])
axs[1].set_xlim([x_min, x_max])
axs[1].plot(t1_MF, Ld, marker='+', markersize=15, markeredgewidth=2, color='red', zorder=500, alpha=.95, label="Critical point", linestyle='none')
axs[1].text(t1_MF+.015, Ld+.01, "(-1.3420, 0.3)", size=16, color='r')

file_names = ["../data/p_star_model/t1_scanning/NEW_MATRIX_BASED_HAMILTONIAN__PS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-2__t1_fin_0.500000__t1_step_0.001000__t2_3.745242__t3_-4.548785__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"]

t2 = 2 * 2*n_nodes/(n_nodes-1) # t2_MF*(CONVERSION_FACTOR)
t3 = 1 * 6*n_nodes**2/(n_nodes**2-3*n_nodes+2) # t3_MF*(CONVERSION_FACTOR)

t1_PS = []
L_avrg_PS = []
L_RMS_PS = []
for file_name in file_names:
    L_vs_t1 = np.genfromtxt(file_name, delimiter=',')

    # Finding averages
    L_avrg = []
    L_RMS = []
    t1 = []
    i = 0
    while i < L_vs_t1.shape[0]:
        prev_t1 = L_vs_t1[i][0]
        L = []
        while i < L_vs_t1.shape[0] and prev_t1 == L_vs_t1[i][0]:
            L.append(L_vs_t1[i][1])
            i+=1
        if np.sqrt(np.array(L).var()) < .2: # Filtering out transition points
            L_avrg.append(np.array(L).mean())
            L_RMS.append(np.sqrt(np.array(L).var()))
            t1.append(prev_t1)
    t1_PS.extend(t1)
    L_avrg_PS.extend(L_avrg)
    L_RMS_PS.extend(L_RMS)

axs[1].plot(t1_PS, L_avrg_PS, color="blue", marker='o', markersize=6, ls='none', label="p-star sim", alpha=.5, zorder=51)

file_names = ["../data/mean_field_model/t1_scanning/NEW_METROPOLIS_MF__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_-2__t1_fin_0.500000__t1_step_0.001000__t2_3.745242__t3_-4.548785__N_DYNAMIC_PAIRS_499500__N_AVRGNG_10__ENSEMBLE_AVRG_0.csv"]

t1_MF = []
L_avrg_MF = []
L_RMS_MF = []
for file_name in file_names:
    L_vs_t1 = np.genfromtxt(file_name, delimiter=',')

    # Finding averages
    L_avrg = []
    L_RMS = []
    t1 = []
    i = 0
    while i < L_vs_t1.shape[0]:
        prev_t1 = L_vs_t1[i][0]
        L = []
        while i < L_vs_t1.shape[0] and prev_t1 == L_vs_t1[i][0]:
            L.append(L_vs_t1[i][1])
            i+=1
        if np.sqrt(np.array(L).var()) < .2: # Filtering out transition points
            L_avrg.append(np.array(L).mean())
            L_RMS.append(np.sqrt(np.array(L).var()))
            t1.append(prev_t1)
    t1_MF.extend(t1)
    L_avrg_MF.extend(L_avrg)
    L_RMS_MF.extend(L_RMS)


axs[1].plot(t1_MF, L_avrg_MF, color="orange", marker='s', markersize=6, ls='none', label="MF sim", alpha=1, zorder=50)

axs[1].legend(loc=[.01, .625], prop={'size': 21})

plt.sca(axs[1])
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)


###########################################################################################

plt.tight_layout()
plt.savefig("../figures/crit_FE_SC.png", dpi=300)

plt.show()
