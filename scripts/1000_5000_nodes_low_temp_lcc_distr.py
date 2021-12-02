###################################################
# This script plots distributions of local        #
# clustering coefficient of the corresponding     #
# 3-star and MF3 models with 1000 and 5000 nodes. #
###################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2,1, figsize=(10., 1.5*6.))

N_nodes = 1000
file_name = "/Users/oleg//sync/study/Northumbria/papers_repos/p_star_and_MF_ERGMs_/data/p_star_model/clustering_distributions/PS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_100000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(float(line.split('\n')[0]))
    
N_lines = np.array(lines).shape[0]
print(np.average(deg_seq_PS))

file_name = "/Users/oleg//sync/study/Northumbria/papers_repos/p_star_and_MF_ERGMs_/data/mean_field_model/clustering_distributions/MF__N_NODES_1000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_100000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_MF = []
for lc in range(0, N_lines):
    deg_seq_MF.append(float(lines[lc].split('\n')[0]))

print(np.average(deg_seq_MF))

#### Finding MF theory prediction for <C> ####
from scipy.optimize import fsolve
Ld = .3
t2_MF = -7500
t3_MF = 4500
t1_MF = 1/2*np.log(Ld/(1-Ld)) - 2*t2_MF*Ld - 3*t3_MF*Ld*Ld

SC = lambda L : L - .5*(1+np.tanh(t1_MF + 2*t2_MF*L + 3*t3_MF*L**2))

L_initial_guess = 0.4
L_sol = fsolve(SC, L_initial_guess)[0]
cc_theor = L_sol*(1-(1-L_sol)**(N_nodes-1)-L_sol*(N_nodes-1)*(1-L_sol)**(N_nodes-2))
print("C_theor=", cc_theor)

# axs[0].set_xlabel(r'$c$', fontsize=26, labelpad=-1)
axs[0].set_ylabel(r'$\langle p(c|\mathbf{A}) \rangle$', fontsize=26)
axs[0].axvline(cc_theor, color='r', label="MF theor", ls='--')
axs[0].hist(deg_seq_PS, 100, color='b', density=True, label="PS sim", histtype='step', linewidth=2)
axs[0].hist(deg_seq_MF, 100, color='orange', density=True, label="MF sim", histtype='step', linewidth=2)
axs[0].legend(fontsize=20, loc="upper right")

plt.sca(axs[0])
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)

##########################################################################################################
N_nodes = 5000

file_name = "/Users/oleg//sync/study/Northumbria/MF_and_p_star_ERGMs/models/data/p_star_model/clustering_distributions/PS__N_NODES_5000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(float(line.split('\n')[0]))
    
N_lines = np.array(lines).shape[0]
print(np.average(deg_seq_PS))


file_name = "/Users/oleg//sync/study/Northumbria/papers_repos/p_star_and_MF_ERGMs_/data/mean_field_model/clustering_distributions/MF__N_NODES_5000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_MF = []
for lc in range(0, N_lines):
    deg_seq_MF.append(float(lines[lc].split('\n')[0]))

print(np.average(deg_seq_MF))

#### Finding MF theory prediction for <C> ####
from scipy.optimize import fsolve
Ld = .3
t2_MF = -7500
t3_MF = 4500
t1_MF = 1/2*np.log(Ld/(1-Ld)) - 2*t2_MF*Ld - 3*t3_MF*Ld*Ld

SC = lambda L : L - .5*(1+np.tanh(t1_MF + 2*t2_MF*L + 3*t3_MF*L**2))

L_initial_guess = 0.4
L_sol = fsolve(SC, L_initial_guess)[0]
cc_theor = L_sol*(1-(1-L_sol)**(N_nodes-1)-L_sol*(N_nodes-1)*(1-L_sol)**(N_nodes-2))
print("C_theor=", cc_theor)

axs[1].set_xlabel(r'$c$', fontsize=26, labelpad=-1)
axs[1].set_ylabel(r'$\langle p(c|\mathbf{A}) \rangle$', fontsize=26, labelpad=-1)
axs[1].axvline(cc_theor, color='r', alpha=1, label="MF theor", ls='--')
axs[1].hist(deg_seq_PS, 100, color='b', density=True, label="PS sim", histtype='step', linewidth=2)
axs[1].hist(deg_seq_MF, 100, color='orange', density=True, label="MF sim", histtype='step', linewidth=2)
axs[1].legend(fontsize=20, loc="upper right")

plt.sca(axs[1])
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)

plt.tight_layout()
plt.savefig('../figures/1000_5000_nodes_low_temp_lcc_distr.png', dpi=300)
plt.show()

fig.clear()
