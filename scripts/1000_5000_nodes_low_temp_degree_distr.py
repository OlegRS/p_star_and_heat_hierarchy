#############################################
# This script plots degree distributions of #
# the corresponding 3-star and MF3 models   #
# with 1000 and 5000 nodes.                 #
#############################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2,1, figsize=(10., 1.5*6.))

N_nodes = 1000
file_name = "../data/p_star_model/degree_distributions/PS__N_NODES_1000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_100000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(int(line.split('\n')[0]))
N_lines = np.array(lines).shape[0]
L_PS = np.average(deg_seq_PS)
print(L_PS)


file_name = "../data/mean_field_model/degree_distributions/MF__N_NODES_1000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_100000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_MF = []
for lc in range(0, N_lines):
    deg_seq_MF.append(int(lines[lc].split('\n')[0]))

L_MF = np.average(deg_seq_MF)
print(L_MF)

# Theoretical
L_MF /= (N_nodes-1)
sigma = np.sqrt((N_nodes-1)*L_MF*(1-L_MF))
x_min, x_max = L_MF*(N_nodes-1) - 3.5*sigma, L_MF*(N_nodes-1) + 3.5*sigma
degrees = np.arange(np.rint(x_min), np.rint(x_max)+1)
p = []
p = np.sqrt(1/(2*np.pi*(N_nodes-1)*L_MF*(1-L_MF)))*np.exp(-(degrees-L_MF*(N_nodes-1))**2/(2*(N_nodes-1)*L_MF*(1-L_MF)))

# axs[0].set_xlabel(r'$k$', fontsize=26)
axs[0].set_ylabel(r'$\langle p(k|\mathbf{A}) \rangle$', fontsize=26)
axs[0].plot(degrees, p, color='r', label="MF theor", ls='--')
axs[0].hist(deg_seq_PS, bins = degrees, color='b', density=True, label="PS sim", histtype='step', linewidth=2)
axs[0].hist(deg_seq_MF, bins = degrees, color='orange', density=True, label="MF sim", histtype='step', linewidth=2)
axs[0].legend(fontsize=20, loc="upper right")
# axs[0].set_title(r"$n=1000;\ t_1=3284.576351,\ t_2=-7500,\ t_3=4500$", fontsize=17)

# axs[0].grid()
axs[0].set_xlim([x_min, x_max])

plt.sca(axs[0])
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)

###########################################################################
N_nodes = 5000
file_name = "../data/p_star_model/degree_distributions/PS__N_NODES_5000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_PS = []
for line in lines:
    deg_seq_PS.append(int(line.split('\n')[0]))
N_lines = np.array(lines).shape[0]
L_PS = np.average(deg_seq_PS)
print(L_PS)


file_name = "../data/mean_field_model/degree_distributions/MF__N_NODES_5000__N_ITERS_PER_LINK_10__t1_MF_3284.576351__t2_MF_-7500__t3_MF_4500__N_AVRGNG_1000__INITIALIZE_RANDOMLY_0.csv"

degs_file = open(file_name, 'r')
lines = degs_file.readlines()
deg_seq_MF = []
for lc in range(0, N_lines):
    deg_seq_MF.append(int(lines[lc].split('\n')[0]))

L_MF = np.average(deg_seq_MF)
print(L_MF)

# Theoretical
L_MF /= (N_nodes-1)
sigma = np.sqrt((N_nodes-1)*L_MF*(1-L_MF))
x_min, x_max = L_MF*(N_nodes-1) - 3.5*sigma, L_MF*(N_nodes-1) + 3.5*sigma
degrees = np.arange(np.rint(x_min), np.rint(x_max)+1)
p = []
p = np.sqrt(1/(2*np.pi*(N_nodes-1)*L_MF*(1-L_MF)))*np.exp(-(degrees-L_MF*(N_nodes-1))**2/(2*(N_nodes-1)*L_MF*(1-L_MF)))


axs[1].set_xlabel(r'$k$', fontsize=26)
axs[1].set_ylabel(r'$\langle p(k|\mathbf{A}) \rangle$', fontsize=26)
axs[1].plot(degrees, p, color='r', label="MF theor", ls='--')
axs[1].hist(deg_seq_PS, bins = degrees, color='b', density=True, label="PS sim", histtype='step', linewidth=2)
axs[1].hist(deg_seq_MF, bins = degrees, color='orange', density=True, label="MF sim", histtype='step', linewidth=2)
axs[1].legend(fontsize=20, loc="upper right")
# axs[1].set_title(r"$n=5000;\ t_1=3284.576351,\ t_2=-7500,\ t_3=4500$", fontsize=17)

# axs[1].grid()
axs[1].set_xlim([x_min, x_max])

plt.sca(axs[1])
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)


plt.tight_layout()
plt.savefig('../figures/1000_5000_nodes_low_temp_degree_distr.png', dpi=300)
plt.show()

fig.clear()
