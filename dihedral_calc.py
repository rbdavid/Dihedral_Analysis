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

pdb = sys.argv[1]
traj_loc = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])
system = sys.argv[5]

nRes = len(res_list.res)

flush = sys.stdout.flush

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
	print '%s' %(string)
        flush()

def summary(nSteps):
	global system

	sum_file = open('%s.dihedral.summary' %(system),'w')
	sum_file.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
	sum_file.write('To recreate this analysis, run this line:\n')
	for i in range(len(sys.argv)):
		sum_file.write('%s ' %(sys.argv[i]))
	sum_file.write('\n\n')
	sum_file.write('output files are written to:\n')
	sum_file.write('	phi_dihedral.%s.dat\n' %(system))
	sum_file.write('	psi_dihedral.%s.dat\n' %(system))
	sum_file.write('	omega_dihedral.%s.dat\n' %(system))
	sum_file.write('	chi1_dihedral.%s.dat\n' %(system))
	sum_file.write('\nNumber of Steps Analyzed: %d\n' %(nSteps))
	sum_file.close()

# ----------------------------------------
# MAIN PROGRAM:

# LOAD IN MOBILE UNIVERSE
u = MDAnalysis.Universe('%s' %(pdb))
u_protein = u.select_atoms('protein')

# MAKE ATOM SELECTIONS FOR RESIDUES (AND DIHEDRALS) OF INTEREST
u_sel = ['']*nRes
for i in range(nRes):
	res_num = res_list.res[i]-1
	u_sel[i] = []
	u_sel[i].append(u_protein.residues[res_num].phi_selection())
	u_sel[i].append(u_protein.residues[res_num].psi_selection())
	u_sel[i].append(u_protein.residues[res_num].omega_selection())
	u_sel[i].append(u_protein.residues[res_num].chi1_selection())

nSel = len(u_sel[0])

# OPEN OUTPUT FILES
out1 = open('phi_dihedral.%s.dat' %(system),'w')
out2 = open('psi_dihedral.%s.dat' %(system),'w')
out3 = open('omega_dihedral.%s.dat' %(system),'w')
out4 = open('chi1_dihedral.%s.dat' %(system),'w')

# CREATE AN OUTPUT FILE LIST...
out_files = []
out_files.append(out1)
out_files.append(out2)
out_files.append(out3)
out_files.append(out4)

nSteps = 0

# BEGINNING TRAJECTORY ANALYSIS
ffprint('Beginning trajectory analysis')
while start <= end:
	ffprint('Loading trajectory %s' %(start))
	u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))
	nSteps += len(u.trajectory)

	# Loop through trajectory
	for ts in u.trajectory:
		# Collect and output dihedral values for each residue of interest
		for i in range(nRes):
			for j in range(nSel):
				if u_sel[i][j] == None:
					out_files[j].write('%10.6f   ' %(0))
				else:
					temp = u_sel[i][j].dihedral.value()
					out_files[j].write('%10.6f   ' %(temp))

		out1.write('\n')
		out2.write('\n')
		out3.write('\n')
		out4.write('\n')

	start += 1

out1.close()
out2.close()
out3.close()
out4.close()

summary(nSteps)

