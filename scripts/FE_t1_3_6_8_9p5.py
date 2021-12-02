###################################################
# This script plots the free energy as a function #
# connectance for different values of t1.         #
###################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

fig, axs = plt.subplots(2, 2, figsize=(10., 6.))

####################################################################################################
t2_MF = -15
t3_MF = 9
t1_MF = 3

# title = "$t_1=$" + str(round(t1_MF, 4)) + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)

##################### Free Energy #######################
L = np.arange(1e-5,1, 1e-5)

F = -2*t1_MF*L - 2*t2_MF*L**2 - 2*t3_MF*L**3 + L*np.log(L/(1-L)) + np.log(1-L)

axs[0,0].plot(L, F, label="Free energy", linewidth=3)

axs[0,0].axvline(0, color='black')
axs[0,0].axvline(1, color='black')
axs[0,0].axhline(0, color='black')

axs[0,0].set_xlim([-.015, 1.015])
axs[0,0].legend(loc='best', prop={'size': 20})

# axs[0,0].set_title(title, fontsize=22)
axs[0,0].set_ylabel(r'$\Phi(\langle L\rangle |{\bf t})$', fontsize=26, labelpad=-.1)

def self_consistency(Ld): # SC equation to find connectance at given t
    return Ld - 1/2*(1 + np.tanh(t1_MF + 2*t2_MF*Ld + 3*t3_MF*Ld**2))

Ld = fsolve(self_consistency, .5)[0]
axs[0,0].plot(Ld, -2*t1_MF*Ld - 2*t2_MF*Ld**2 - 2*t3_MF*Ld**3 + Ld*np.log(Ld/(1-Ld)) + np.log(1-Ld),  marker='o', fillstyle='none', markersize=12, color='black', zorder=500, alpha=1, markeredgewidth=4)

#########################################################################################################
t2_MF = -15
t3_MF = 9
t1_MF = 6

# title = "$t_1=$" + str(round(t1_MF, 4)) + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)

##################### Free Energy #######################
L = np.arange(1e-5,1, 1e-5)

F = -2*t1_MF*L - 2*t2_MF*L**2 - 2*t3_MF*L**3 + L*np.log(L/(1-L)) + np.log(1-L)

axs[0,1].plot(L, F, label="Free energy", linewidth=3)

axs[0,1].axvline(0, color='black')
axs[0,1].axvline(1, color='black')
axs[0,1].axhline(0, color='black')

axs[0,1].set_xlim([-.015, 1.015])
# axs[0,1].legend(loc='best', prop={'size': 20})

# axs[0,1].set_title(title, fontsize=22)
# axs[0,1].set_ylabel(r'$\Phi(\langle L\rangle |{\bf t})$', fontsize=22)

def self_consistency(Ld): # SC equation to find connectance at given t
    return Ld - 1/2*(1 + np.tanh(t1_MF + 2*t2_MF*Ld + 3*t3_MF*Ld**2))

Ld = fsolve(self_consistency, .3)[0] # Lower solution
axs[0,1].plot(Ld, -2*t1_MF*Ld - 2*t2_MF*Ld**2 - 2*t3_MF*Ld**3 + Ld*np.log(Ld/(1-Ld)) + np.log(1-Ld), marker='^', fillstyle='none', markersize=12, color='forestgreen', zorder=500, alpha=1, markeredgewidth=3)

Lu = fsolve(self_consistency, .99)[0] # Upper solution
axs[0,1].plot(Lu, -2*t1_MF*Lu - 2*t2_MF*Lu**2 - 2*t3_MF*Lu**3 + Lu*np.log(Lu/(1-Lu)) + np.log(1-Lu), marker='v', fillstyle='none', markersize=12, color='brown', zorder=500, alpha=1, markeredgewidth=3)
########################################################################################
t2_MF = -15
t3_MF = 9
t1_MF = 8

# title = "$t_1=$" + str(round(t1_MF, 4)) + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)

##################### Free Energy #######################
L = np.arange(1e-5,1, 1e-5)

F = -2*t1_MF*L - 2*t2_MF*L**2 - 2*t3_MF*L**3 + L*np.log(L/(1-L)) + np.log(1-L)

axs[1,0].plot(L, F, label="Free energy", linewidth=3)

axs[1,0].axvline(0, color='black')
axs[1,0].axvline(1, color='black')
axs[1,0].axhline(0, color='black')

axs[1,0].set_xlim([-.015, 1.015])
# axs[1,0].legend(loc='best', prop={'size': 20})

# axs[1,0].set_title(title, fontsize=22)
axs[1,0].set_xlabel(r'$\langle L\rangle$', fontsize=26, labelpad=-.1)
axs[1,0].set_ylabel(r'$\Phi(\langle L\rangle |{\bf t})$', fontsize=26, labelpad=-.1)

def self_consistency(Ld): # SC equation to find connectance at given t
    return Ld - 1/2*(1 + np.tanh(t1_MF + 2*t2_MF*Ld + 3*t3_MF*Ld**2))

Ld = fsolve(self_consistency, .3)[0] # Lower solution
axs[1,0].plot(Ld, -2*t1_MF*Ld - 2*t2_MF*Ld**2 - 2*t3_MF*Ld**3 + Ld*np.log(Ld/(1-Ld)) + np.log(1-Ld), marker='<', fillstyle='none', markersize=12, color='brown', zorder=500, alpha=1, markeredgewidth=3)

Lu = fsolve(self_consistency, .99)[0] # Upper solution
axs[1,0].plot(Lu, -2*t1_MF*Lu - 2*t2_MF*Lu**2 - 2*t3_MF*Lu**3 + Lu*np.log(Lu/(1-Lu)) + np.log(1-Lu), marker='>', fillstyle='none', markersize=12, color='forestgreen', zorder=500, alpha=1, markeredgewidth=3)
#########################################################################################################
t2_MF = -15
t3_MF = 9
t1_MF = 9.5

# title = "$t_1=$" + str(round(t1_MF, 4)) + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)

##################### Free Energy #######################
L = np.arange(1e-5,1, 1e-5)

F = -2*t1_MF*L - 2*t2_MF*L**2 - 2*t3_MF*L**3 + L*np.log(L/(1-L)) + np.log(1-L)

axs[1,1].plot(L, F, label="Free energy", linewidth=3)

axs[1,1].axvline(0, color='black')
axs[1,1].axvline(1, color='black')
axs[1,1].axhline(0, color='black')

axs[1,1].set_xlim([-.015, 1.015])
# axs[1,1].legend(loc='best', prop={'size': 20})

# axs[1,1].set_title(title, fontsize=22)
axs[1,1].set_xlabel(r'$\langle L\rangle$', fontsize=26, labelpad=-.1)
# axs[1,1].set_ylabel(r'$\Phi(\langle L\rangle |{\bf t})$', fontsize=22)

def self_consistency(Ld): # SC equation to find connectance at given t
    return Ld - 1/2*(1 + np.tanh(t1_MF + 2*t2_MF*Ld + 3*t3_MF*Ld**2))

Ld = fsolve(self_consistency, .5)[0]
axs[1,1].plot(Ld, -2*t1_MF*Ld - 2*t2_MF*Ld**2 - 2*t3_MF*Ld**3 + Ld*np.log(Ld/(1-Ld)) + np.log(1-Ld),  marker='*', fillstyle='none', markersize=14, color='black', zorder=500, alpha=1, markeredgewidth=3)
###############################################################################################################

axs[0,0].set_yticks([-1, 0, 1, 2, 3, 4, 5, 6])
axs[0,0].set_yticklabels([-1, 0, 1, 2, 3, 4, 5, 6], fontsize=17)
axs[0,1].set_yticks([-2, -1.5, -1, -.5, 0])
axs[0,1].set_yticklabels([-2, -1.5, -1, -.5, 0], fontsize=17)
axs[1,0].set_yticks([-4, -3, -2, -1, 0])
axs[1,0].set_yticklabels([-4, -3, -2, -1, 0], fontsize=17)
axs[1,1].set_yticks([-7, -6, -5, -4, -3, -2, -1, 0])
axs[1,1].set_yticklabels([-7, -6, -5, -4, -3, -2, -1, 0], fontsize=17)


axs[0,0].set_xticks([0, .2, .4, .6, .8, 1])
axs[0,0].set_xticklabels([])
axs[0,1].set_xticks([0, .2, .4, .6, .8, 1])
axs[0,1].set_xticklabels([])
axs[1,0].set_xticks([0, .2, .4, .6, .8, 1])
axs[1,0].set_xticklabels([0, .2, .4, .6, .8, 1], fontsize=17)
axs[1,1].set_xticks([0, .2, .4, .6, .8, 1])
axs[1,1].set_xticklabels([0, .2, .4, .6, .8, 1], fontsize=17)

plt.tight_layout()
plt.savefig("../figures/FE_t1_3_6_8_9p5.png", dpi=300)

plt.show()
