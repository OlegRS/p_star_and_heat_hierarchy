######################################################
# This script plots free energy of 3-star model as a #
# function of connectance for the given t1, t2, t3.  #
######################################################
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import fsolve


n_nodes = 1000
n_pairs = n_nodes*(n_nodes-1)/2

fig, ax = plt.subplots(figsize=(10., 6.))

t2_MF = 2
t3_MF = 1
t1_MF = -3

title = "$t_1=$" + str(t1_MF) + "$;\ t_2=$" + str(t2_MF) + "$;\ t_3=$" + str(t3_MF)

dL = 1e-7
L = np.arange(dL, 1+2*dL, dL)

F = -2*t1_MF*L - 2*t2_MF*L**2 - 2*t3_MF*L**3 + L*np.log(L/(1-L)) + np.log(1-L)
plt.plot(L, F, label="Free energy")

plt.axvline(0, color='black')
plt.axhline(0, color='black')

plt.xlim([0, 1])

# plt.title(title, fontsize=22)
plt.xlabel(r'$\langle L\rangle$', fontsize=26, labelpad=-8)
plt.ylabel(r'$\Phi(\langle L\rangle | {\bf t})$', fontsize=26, labelpad=-1)

# plt.grid()

plt.legend(loc='best', prop={'size': 22})

def self_consistency(Ld): # SC equation to find connectance at given t
    return Ld - 1/2*(1 + np.tanh(t1_MF + 2*t2_MF*Ld + 3*t3_MF*Ld**2))

Ld_0 = fsolve(self_consistency, 0)[0]
Ld_1 = fsolve(self_consistency, .99)[0]
print('Ld_0=', Ld_0)
print('Ld_1=', Ld_1)

plt.plot(Ld_0, -2*t1_MF*Ld_0 - 2*t2_MF*Ld_0**2 - 2*t3_MF*Ld_0**3 + Ld_0*np.log(Ld_0/(1-Ld_0)) + np.log(1-Ld_0), marker='+', markersize=15, color='red', zorder=500, markeredgewidth=2)
plt.plot(Ld_1, -2*t1_MF*Ld_1 - 2*t2_MF*Ld_1**2 - 2*t3_MF*Ld_1**3 + Ld_1*np.log(Ld_1/(1-Ld_1)) + np.log(1-Ld_1), marker='+', markersize=15, color='red', zorder=500, markeredgewidth=2)

plt.yticks(fontsize=17)
plt.xticks(fontsize=17)

##### ZOOM ######
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
axins = zoomed_inset_axes(ax, 29, bbox_to_anchor=(.4, .47+.05, .6, .52+.05), bbox_transform=ax.transAxes, loc="lower left")
axins.plot(L, F)
axins.axhline(0, color='black')
# axins.grid()
x1, x2, y1, y2 = 0, .01, -0.005, 0.005 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits

axins.plot(Ld_0, -2*t1_MF*Ld_0 - 2*t2_MF*Ld_0**2 - 2*t3_MF*Ld_0**3 + Ld_0*np.log(Ld_0/(1-Ld_0)) + np.log(1-Ld_0), marker='+', markersize=15, color='red', zorder=500, markeredgewidth=2)

axins.tick_params(labelsize=14)
axins.tick_params(axis='both', which='minor', labelsize=14)
axins.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
axins.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
axins.yaxis.offsetText.set_fontsize(14)
axins.xaxis.offsetText.set_fontsize(14)

from mpl_toolkits.axes_grid1.inset_locator import mark_inset
mark_inset(ax, axins, loc1=2, loc2=3, fc="none", ec="0.5")

axins = zoomed_inset_axes(ax, 290, bbox_to_anchor=(.435, .095+.01, .6, .5+.01), bbox_transform=ax.transAxes, loc="lower left") # zoom-factor: 5
axins.plot(L, F)
axins.axhline(0, color='black')
# axins.grid()

x1, x2, y1, y2 = .9988, 1, -.0005, .0005 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits

axins.tick_params(labelsize=14)
axins.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
axins.ticklabel_format(axis='x', style='plain', scilimits=(0,0))
axins.yaxis.offsetText.set_fontsize(14)
axins.xaxis.offsetText.set_fontsize(14)
axins.set_xticks(np.arange(0.999, 1, step=0.0005))

plt.plot(Ld_1, -2*t1_MF*Ld_1 - 2*t2_MF*Ld_1**2 - 2*t3_MF*Ld_1**3 + Ld_1*np.log(Ld_1/(1-Ld_1)) + np.log(1-Ld_1), marker='+', markersize=15, color='red', zorder=500, markeredgewidth=2)

mark_inset(ax, axins, loc1=1, loc2=4, fc="none", ec="0.5")

plt.tight_layout()
plt.savefig("../figures/free_energy.png", dpi=300)

plt.show()
