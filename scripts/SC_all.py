########################################################
# This script plots the MF self-consistency equation   #
# together with the results of Monte Carlo simulations #
# for 3-star and MF3 models. The marked points show    #
# the positions of the free energy minima from         #
# ../figures/FE_t1_3_6_8_9p5.png                       #
########################################################

import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

fig, ax = plt.subplots(1, 1, figsize=(10., 6.))
t2_MF = -15
t3_MF = 9
t1_MF = 3
title = r"$n=10^3$" + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)
n_nodes = 1000
################## t1 DEPENDENCE ########################
x_min, x_max =2, 10
delta_t1 = 0.05
delta_L = 0.00001
t1_range = np.arange(x_min, x_max+delta_t1, delta_t1)
L_range = np.arange(0, 1+2*delta_L, delta_L)
L, T1 = np.meshgrid(L_range, t1_range)
SC = L - 1/2*(1 + np.tanh(T1 +2*t2_MF*L + 3*t3_MF*L**2)) #Self-consistency
CS = ax.contour(T1, L, SC, [0], colors='cyan', zorder=100, linewidths=1.5)
CS.collections[0].set_label('MF solution')

def self_consistency(Ld): # SC equation to find connectance at given t
    return Ld - 1/2*(1 + np.tanh(t1_MF + 2*t2_MF*Ld + 3*t3_MF*Ld**2))
Ld = fsolve(self_consistency, .5)[0]
ax.plot(t1_MF, Ld, marker='o', fillstyle='none', markersize=12, color='black', zorder=500, alpha=1, markeredgewidth=4)

t1_MF = 6
Ld = fsolve(self_consistency, .5)[0]
ax.plot(t1_MF, Ld, marker='^', fillstyle='none', markersize=12, color='forestgreen', zorder=500, alpha=1, markeredgewidth=3)
Ld = fsolve(self_consistency, .99)[0]
ax.plot(t1_MF, Ld, marker='v', fillstyle='none', markersize=12, color='brown', zorder=500, alpha=1, markeredgewidth=3)

t1_MF = 8
Ld = fsolve(self_consistency, .5)[0]
ax.plot(t1_MF, Ld, marker='<', fillstyle='none', markersize=12, color='brown', zorder=500, alpha=1, markeredgewidth=3)
Ld = fsolve(self_consistency, .99)[0]
ax.plot(t1_MF, Ld, marker='>', fillstyle='none', markersize=12, color='forestgreen', zorder=500, alpha=1, markeredgewidth=3)

t1_MF = 9.5
Ld = fsolve(self_consistency, .99)[0]
ax.plot(t1_MF, Ld, marker='*', fillstyle='none', markersize=14, color='black', zorder=500, alpha=1, markeredgewidth=3)


ax.set_ylabel(r'$\langle L\rangle$', fontsize=26, labelpad=-6)
ax.set_xlabel(r'$t_1$', fontsize=26, labelpad=-3)
ax.set_xlim([x_min, x_max])
ax.set_ylim([-.015, 1.025])
ax.axvline(0, color="black")
ax.axhline(0, color="black")
# ax.grid()

file_names = ["../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_2.000000__t1_fin_10.100000__t1_step_0.010000__t2_MF_-15.000000__t3_MF_9.000000__N_AVRGNG_10__RANDOM_INITIALISATION_0.csv",
              "../data/p_star_model/t1_scanning/PS__CONNECTANCE_and_STARS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_10.000000__t1_fin_1.990000__t1_step_-0.010000__t2_MF_-15.000000__t3_MF_9.000000__N_AVRGNG_10__RANDOM_INITIALISATION_0.csv"
]

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

ind = t1_PS.index(8.48) # Removing transition point
del t1_PS[ind]
del L_avrg_PS[ind]
del L_RMS_PS[ind]

# ax.errorbar(t1_PS, L_avrg_PS, xerr=None, yerr=L_RMS_PS, color="blue", marker='o', markersize=4, ls='none', label="PS simulation", alpha=.5, elinewidth=2, capsize=2, zorder=51)
ax.plot(t1_PS, L_avrg_PS, color="blue", marker='o', markersize=4, ls='none', label="PS simulation", alpha=.5, zorder=51)

file_names = ["../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_2.000000__t1_fin_10.100000__t1_step_0.010000__t2_MF_-15.000000__t3_MF_9.000000__N_DYNAMIC_PAIRS_1__N_AVRGNG_10__RANDOM_INITIALISATION_0.csv",
              "../data/mean_field_model/t1_scanning/MF__CONNECTANCE_and_STARS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_start_10.000000__t1_fin_1.990000__t1_step_-0.010000__t2_MF_-15.000000__t3_MF_9.000000__N_DYNAMIC_PAIRS_1__N_AVRGNG_10__RANDOM_INITIALISATION_0.csv"
]

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

# ax.errorbar(t1_MF, L_avrg_MF, xerr=None, yerr=L_RMS_MF, color="orange", marker='s', markersize=4, ls='none', label="MF simulation", alpha=1, elinewidth=2, capsize=2, zorder=50)
ax.plot(t1_MF, L_avrg_MF, color="orange", marker='s', markersize=6, ls='none', label="MF simulation", alpha=1, zorder=50)

#### Thermodynamic solution ####
x_min_ = 7
x_max_ = 8

T1_MF = np.arange(x_min_, x_max_, 0.01)
L_MIN = []
L_EXPECTED = []
T1_MF_ = []
for t1_MF in T1_MF:
    def SC_equation(L):
        return L - 1/2*(1 + np.tanh(t1_MF + 2*t2_MF*L + 3*t3_MF*L**2))
    L_ = np.array([fsolve(SC_equation, 0)[0], fsolve(SC_equation, .8)[0], fsolve(SC_equation, .3)[0], fsolve(SC_equation, .999)[0], fsolve(SC_equation, 1.001)[0], fsolve(SC_equation, .99)[0], fsolve(SC_equation, 1)[0]])
    def F(L):
        return -2*(t1_MF*L + t2_MF*L**2 + t3_MF*L**3) + L*np.log(L/(1-L)) + np.log(1-L)
    L0 = L_[~(np.triu(np.abs(L_[:,None] - L_) < 0.001,1)).any(0)] # Distinct solution of SC equation
    L = np.array([l for l in L0 if (l<1 and l>0)])        
    L_min = L_[list(F(L_)).index(min(F(L_)))]
    L_MIN.append(L_min)

ax.plot(T1_MF, L_MIN, color="red", label='Thermodynamic solution', ls='--')
handles, labels = ax.get_legend_handles_labels()
ax.legend([handles[2], handles[3], handles[1], handles[0]], [labels[2],labels[3],labels[1],labels[0]], prop={'size': 20}, loc=[.01, .5])
# ax.set_title(title, fontsize=22)

plt.yticks(fontsize=17)
plt.xticks(fontsize=17)

plt.tight_layout()
plt.savefig("../figures/SC_all.png", dpi=300)

plt.show()
