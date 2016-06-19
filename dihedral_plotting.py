#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:


# ----------------------------------------
# PREAMBLE:

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import res_list

nullfmt = NullFormatter()

phi_file = sys.argv[1]
psi_file = sys.argv[2]
omega_file = sys.argv[3]
system = sys.argv[4]

nRes = len(res_list.res)

phi_data = np.loadtxt(phi_file)
psi_data = np.loadtxt(psi_file)
omega_data = np.loadtxt(omega_file)

my_cmap = plt.cm.get_cmap('jet')
my_cmap.set_under('w')

for i in range(nRes):
# IF ONLY ONE RESIDUE IN RES_LIST FILE
#	counts, xedges, yedges, image = plt.hist2d(phi_data[:], psi_data[:], bins=100, cmap=my_cmap,vmin = 0.1)
# IF MULTIPLE RESIDUES IN RES_LIST FILE
	counts, xedges, yedges, image = plt.hist2d(phi_data[:,i], psi_data[:,i], bins=100, cmap=my_cmap,vmin = 0.1)
	cb1 = plt.colorbar()
	cb1.set_label('Frequency')
	plt.xlabel('Phi Dihedral')
	plt.ylabel('Psi Dihedral')
	plt.xlim((-180,180))
	plt.ylim((-180,180))
	plt.xticks(np.arange(-180,181,90))
	plt.yticks(np.arange(-180,181,90))
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.savefig('%03d.%s.rama_plot.png' %(res_list.res[i],system))
	plt.close()
	counts = []
	xedges = []
	yedges = []
	image = []

# IF ONLY ONE RESIDUE IN RES_LIST FILE
#	plt.plot(omega_data[:],'k.')
# IF MULTIPLE RESIDUES IN RES_LIST FILE
	plt.plot(omega_data[:,i],'k.')
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel(r'Time (timesteps)', size=12)
	plt.ylabel(r'Omega Dihedral', size=12)
	plt.ylim((-180,180))
	plt.yticks(np.arange(-180,181,90))
	plt.savefig('%03d.%s.omega_dihedral.png' %(res_list.res[i],system))
	plt.close()

