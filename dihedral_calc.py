#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:


# ----------------------------------------
# PREAMBLE:

import numpy as np
import sys
import os
import MDAnalysis
from MDAnalysis.analysis.align import *
import res_list
from distance_functions import *

ref_frame = sys.argv[1]
pdb = sys.argv[2]
traj_loc = sys.argv[3]
start = int(sys.argv[4])
end = int(sys.argv[5])
system = sys.argv[6]

alignment = 'protein and name CA and (resid 20:25 or resid 50:55 or resid 73:75 or resid 90:94 or resid 112:116 or resid 142:147 or resid 165:169 or resid 190:194 or resid 214:218 or resid 236:240 or resid 253:258 or resid 303:307)'

nRes = len(res_list.res)

flush = sys.stdout.flush

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
	print '%s' %(string)
        flush()

# ----------------------------------------
# MAIN PROGRAM:
# LOAD IN REFERENCE STRUCTURE 
ref = MDAnalysis.Universe('%s' %(ref_pdb))
ref_align = ref.select_atoms(alignment)
ref_all = ref.select_atoms('all')
ref_backbone = ref.select_atoms('backbone')
ref_all.translate(-ref_backbone.center_of_mass())
pos0 = ref_align.coordinates()

# LOAD IN MOBILE UNIVERSE
u = MDAnalysis.Universe('%s' %(pdb))
u_align = u.select_atoms(alignment)
u_all = u.select_atoms('all')
u_backbone = u.select_atoms('backbone')

if len(u_align) != len(ref_align):
	ffprint('Alignment atom selections do not have the same number of atoms.')
	sys.exit()

# MAKE ATOM SELECTIONS FOR BOTH MOBILE UNIVERSES
u_sel = ['']*nRes
for i in range(nRes):
	res_num = res_list.res[i]-1
	u_sel[i] = []
	u_sel[i].append(u_all.residues[res_num].phi_selection())
	u_sel[i].append(u_all.residues[res_num].psi_selection())

# OPEN OUTPUT FILES
out1 = open('phi_dihedral.dat')
out2 = open('psi_dihedral.dat')

nSteps = 0

# BEGINNING TRAJECTORY ANALYSIS
ffprint('Beginning trajectory analysis')
while start <= end:
	ffprint('Loading trajectory %s' %(start))
	u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))
	nSteps += len(u.trajectory)

	# Loop through trajectory
	for ts in u.trajectory:
		# Align to reference (moves COM of backbone to origin)
		u_all.translate(-u_backbone.center_of_mass())
		# Calculate the rotational matrix to align u to the ref
		R, rmsd = rotation_matrix(u_align.coordinates(), pos0)
		# Apply rotation matrix to atoms within u
		u_all.rotate(R)

		for i in range(nRes):
			phi = u_sel[i][0].dihedral.value()
			psi = u_sel[i][1].dihedral.value()

			out1.write('%10.6f   ' %(phi))
			out2.write('%10.6f   ' %(psi))

		out1.write('\n')
		out2.write('\n')

	start += 1

out1.close()
out2.close()

print 'Analyzed %d steps' %(nSteps)

