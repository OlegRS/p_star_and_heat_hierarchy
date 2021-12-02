##################################################
# This script plots the degree distributions     #
# together with the distributions of local       #
# clustering coefficients for the corresponding  #
# 3-star and MF models.                          #
##################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as ss

fig, axs = plt.subplots(2, 2, figsize=(10., 6.))

N_nodes = 2000
file_name = "../data/p_star_model/degree_distributions/PS__N_NODES_2000__N_ITERS_PER_LINK_10__t1_-0.500000__t2_0.600000__t3_0.600000__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(int(line.split('\n')[0]))
N_lines = np.array(lines).shape[0]
L_PS = np.average(deg_seq_PS)
print(L_PS)


file_name = "../data/mean_field_model/degree_distributions/MF__N_NODES_2000__N_ITERS_PER_LINK_10__t1_-0.500000__t2_0.600000__t3_0.600000__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_MF = []
for lc in range(0, N_lines):
    deg_seq_MF.append(int(lines[lc].split('\n')[0]))

L_MF = np.average(deg_seq_MF)
print(L_MF)

# axs[0,0].set_title(r"$t_1=-0.5,\ t_2=0.3,\ t_3=0.1,\ n=2\cdot 10^3$", fontsize=15)
# Theoretical
L_MF /= (N_nodes-1)
sigma = np.sqrt((N_nodes-1)*L_MF*(1-L_MF))
x_min, x_max = L_MF*(N_nodes-1) - 3.5*sigma, L_MF*(N_nodes-1) + 3.5*sigma
degrees = np.arange(np.rint(x_min), np.rint(x_max)+1)
p = []
p = np.sqrt(1/(2*np.pi*(N_nodes-1)*L_MF*(1-L_MF)))*np.exp(-(degrees-L_MF*(N_nodes-1))**2/(2*(N_nodes-1)*L_MF*(1-L_MF)))
axs[0,0].plot(degrees, p, color='r', label="MFT", linestyle='--')
axs[0,0].hist(deg_seq_PS, bins = degrees, color='b', density=True, label="PS", histtype='step', linewidth=2)
axs[0,0].hist(deg_seq_MF, bins = degrees, color='orange', density=True, label="MF", histtype='step', linewidth=2)
axs[0,0].legend(fontsize=13.5, loc=[.67, .58])
axs[0,0].set_yticks([0, .005, .01, .015])
plt.sca(axs[0, 0])
plt.xticks([720, 760, 800, 840], fontsize=17)
plt.yticks(fontsize=17)


# axs[0,0].grid()
axs[0,0].set_xlim([x_min, x_max])

axs[1,0].plot(degrees, p, color='r', linestyle='--')
axs[1,0].hist(deg_seq_PS, bins = degrees, color='b', density=True, histtype='step', linewidth=2)
axs[1,0].hist(deg_seq_MF, bins = degrees, color='orange', density=True, histtype='step', linewidth=2)
axs[1,0].set_xscale('log')
axs[1,0].set_yscale('log')
axs[1,0].minorticks_off()
axs[1,0].set_xticks([720, 760, 800, 840])
axs[1,0].set_xticklabels([720, 760, 800, 840], fontsize=17)
plt.sca(axs[1, 0])
plt.yticks(fontsize=17)

# axs[1,0].grid()
axs[1,0].set_xlim([x_min, x_max])


file_name = "../data/p_star_model/clustering_distributions/PS__N_NODES_2000__N_ITERS_PER_LINK_10__t1_-0.500000__t2_0.600000__t3_0.600000__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(float(line.split('\n')[0]))
    
N_lines = np.array(lines).shape[0]
print(np.average(deg_seq_PS))




file_name = "../data/mean_field_model/clustering_distributions/MF__N_NODES_2000__N_ITERS_PER_LINK_10__t1_-0.500000__t2_0.600000__t3_0.600000__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_MF = []
for lc in range(0, N_lines):
    deg_seq_MF.append(float(lines[lc].split('\n')[0]))

print(np.average(deg_seq_MF))


#### Finding MF theory prediction for <C> ####
from scipy.optimize import fsolve
t1 = -.5
t2 = .3
t3 = .1

SC = lambda L : L - .5*(1+np.tanh(t1 + 2*t2*L + 3*t3*L**2))

L_initial_guess = 0.4
L_sol = fsolve(SC, L_initial_guess)[0]
cc_theor = L_sol*(1-(1-L_sol)**(N_nodes-1)-L_sol*(N_nodes-1)*(1-L_sol)**(N_nodes-2))
print("C_theor=", cc_theor)



p = np.sqrt(1/(2*np.pi*(N_nodes-1)*L_MF*(1-L_MF)))*np.exp(-(degrees-L_MF*(N_nodes-1))**2/(2*(N_nodes-1)*L_MF*(1-L_MF)))
axs[0,1].hist(deg_seq_PS, 100, color='b', density=True, label="PS", histtype='step', linewidth=2)
axs[0,1].hist(deg_seq_MF, 100, color='orange', density=True, label="MF", histtype='step', linewidth=2)
# axs[0,1].grid()
# axs[0,1].set_title(r"$t_1=-0.5,\ t_2=0.3,\ t_3=0.1,\ n=2\cdot 10^3$", fontsize=15)
axs[0,1].axvline(cc_theor, color='r', label="MFT", ls='--')
axs[0,1].legend(fontsize=13.5, loc=[.67, .58])
plt.sca(axs[0, 1])
plt.xticks([.390, .395], fontsize=17)
# plt.xticks([.389, .391, .393, .395], fontsize=17)
plt.yticks(fontsize=17)



axs[1,1].axvline(cc_theor, color='r', ls='--')
axs[1,1].hist(deg_seq_PS, 100, color='b', density=True, histtype='step', linewidth=2)
axs[1,1].hist(deg_seq_MF, 100, color='orange', density=True, histtype='step', linewidth=2)
axs[1,1].set_xscale('log')
axs[1,1].set_yscale('log')
axs[1,1].minorticks_off()
axs[1,1].set_xticks([.390, .395])
plt.sca(axs[1, 1])
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)



# axs[1,1].grid()
axs[1,1].set_xlim([0.389, 0.396])
axs[0,1].set_xlim([0.389, 0.396])


axs[0,0].set_ylabel(r'$\langle p(k|\mathbf{A}) \rangle$', fontsize=26)
axs[1,0].set_ylabel(r'$\langle p(k|\mathbf{A}) \rangle$', fontsize=26)
axs[0,1].set_ylabel(r'$\langle p(c|\mathbf{A}) \rangle$', fontsize=26)
axs[1,1].set_ylabel(r'$\langle p(c|\mathbf{A}) \rangle$', fontsize=26)
axs[1,0].set_xlabel(r'$k$', fontsize=26)
axs[1,1].set_xlabel(r'$c$', fontsize=26)

plt.tight_layout()

plt.savefig('../figures/deg_and_lcc_distrib.png', dpi=300)
plt.show()

fig.clear()
