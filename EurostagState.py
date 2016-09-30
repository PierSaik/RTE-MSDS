# -*- coding: utf-8 -*-

# Copyright (c) 2016, Pierre Saikaly  (saikalypierre@gmail.com)

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this file, 
# You can obtain one at http://mozilla.org/MPL/2.0/.

#===========================#
# created on 07 sept 2016
#===========================#

# Import Python dependencies :
# ----------------------------
import os
import sys
import time
import argparse

# Script :
# --------
def CheckEurostagState(folder):
	"""
		Reads eurostagOutputFile and check if the situation at the beginning is stable
		or unstable 
		Inputs :
			- eurostagOutputFile : file.out containing the informations
		Outputs :
			- 0 -> Situation is stable
			- 1 -> Situation is unstable
		Used in : 
			- main  
	"""

	etat   = "ETAT D'EQUILIBRE A  0.10000D-02 VERIFIE POUR LES EQUATIONS MACHINE"
	etat2  = "ETAT D'EQUILIBRE A  0.10000D-02 NON VERIFIE DANS LES EQUATIONS SUIVANTES"
	eurostagOutputFile = os.path.join(folder,"sim_pre_fault.out")



	# Checking Eurostag status :
	# --------------------------
	if etat in open(eurostagOutputFile, 'U').read():
		# Eurostag returns stable
		status = 0
	else:
		# Eurostag returns unstable :
		# ---------------------------
		equi_value = []
		with open(eurostagOutputFile, 'U') as outfile:
			lines = outfile.readlines()
			for i, line in enumerate(lines):
				if etat2 in line:
					k = i+5
					while (lines[k]!="1\n" and lines[k]!="\n"):
						machine_name = lines[k].split()[0]
						value_line = lines[k].split()[-2]
						value_real = float(value_line.split("D")[0])*10**int(value_line.split("D")[1])
						if (abs(value_real) > 0.01):
							with open('ecart_groupe.csv','a') as egf:
								egf.write(folder + ';' + machine_name + ';' + str(value_real) + '\n')
						equi_value.append(abs(value_real))
						k += 1
		# Checking tolerance :
		# --------------------
		if (len(equi_value) > 0):
			if (max(equi_value) < 0.01):
				status = 0       # Stable state within tolerance
			else:
				status = 1 # Unstable state

		else:
			status = 0       # Stable state

	return status

def main():
	"""
		Get list of eurostag simulation folders and check if thoses situations are
		stable or not. The outputs are stored in eurostag_status.txt file 
	"""

	# Creating output file :
	# ----------------------
	print "Beginning writing outputs for folders in directory..."

	with open('eurostag_status.txt','w') as f:
		f.write("# " + time.strftime("%d/%m/%y") + "\n")

	with open('ecart_groupe.csv','w') as egf:
		egf.write(time.strftime("%d/%m/%y") + "\n")
		egf.write("Situations;groupe;valeur\n")

	# Checking folders and writing output :
	# -------------------------------------
	with open('eurostag_status.txt','a') as f:
		for folder in os.listdir('.'):
			if folder.startswith('itesla_eurostag_stabilization_'):
				f.write(folder + " " + str(CheckEurostagState(folder)) + '\n')
	
	print "Done writing outputs in eurostag_status.txt"
	print "Done analysis in ecart_groupe.csv"

if __name__=='__main__':

	parser = argparse.ArgumentParser(description='Check eurostag status of working directory of specific folder')
	parser.add_argument('-p','--eur_path', help='path of eurostag folder', required=False)
	args = parser.parse_args()

	if (args.eur_path):
		with open('ecart_groupe.csv','w') as egf:
			egf.write(time.strftime("%d/%m/%y") + "\n")
			egf.write("Situations;groupe;valeur\n")

		print args.eur_path, ' ', CheckEurostagState(args.eur_path)
		print "Done analysis in ecart_groupe.csv"
	else:
		main()