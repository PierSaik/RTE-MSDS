# -*- coding: utf-8 -*-

# Copyright (c) 2016, Pierre Saikaly  (saikalypierre@gmail.com)
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#===========================#
# created on 07 july 2016
#===========================#

# Import Python dependencies :
# ----------------------------
import imp
import os
import sys
import glob
import numpy as np
import time

# Script :
# --------

def FfindDict(file_dta):
	"""
		Find the dictionnary that corresponds to the dta file.
		This is necessary to get the reference names of the generators, as 
		the internal name and the reference name might change from file to file
		Inputs : 
			- file_dta : dta file path
		Outputs :
			- file_csv : dictionnary file path
		Used in :
			- main
	"""

	for file_csv in glob.glob('*.csv'):
		if (file_csv.split('_')[3] == file_dta.split('_')[3]):
			print "dta - dict Match found !"
			break
	return file_csv

def FfindRealGenName(MachineName,file_dictPath):
	"""
		From the MachineName and the dictionnary the Reference name is extracted
		Inputs : 
			- MachineName   : This is the name used in the simulation
			- file_dictPath : Path of the dictionnary file
		Output :
			- MachineDictName : Reference name
		Used in :
			- FreadWritedtaFile
	"""

	with open(file_dictPath) as file_dict:
		data_dict = np.loadtxt(file_dict, dtype='string', delimiter=';')

	MachineDictName = data_dict[MachineName==data_dict[:,1],0]

	return MachineDictName

def FreadWritedtaFile(file_dtaPath, file_dict, generators):
	"""
		Append to sim.dta the information about generators from the file_dataPath.
		It also checks if the generators have not already been included by checking 
		the generators list input
		Inputs : 
			- file_dtaPath : dta path of file from which generators are extracted
			- file_dict    : dictionnary path of corresponding dta file
			- generators   : list of generators included already
		Outputs :
			- generators : the list of generators updated with the added generators
		Used in :
			- main
	"""

	with open(file_dtaPath) as file_dta:
		lines = file_dta.readlines()

		#Looking for machines :
		#----------------------
		for i, line in enumerate(lines):
			if line.startswith("M2       U") or line.startswith("M2       S"):

				# Getting machine name and real name :
				# ------------------------------------
				MachineName     = lines[i+1].split(None, 1)[0]
				MachineDictName = FfindRealGenName(MachineName,file_dict)

				# Checking if generator is new :
				# ------------------------------
				if MachineDictName[0] not in generators[:,0]:
					generators = np.vstack((generators,[MachineDictName[0], MachineName]))

					# Adding Machine to sim.dta :
					# ---------------------------
					with open('./rez/sim.dta', 'a') as f:
						#Writing Machine information :
						#-----------------------------
						f.write(lines[i])
						f.write(lines[i+1])
						f.write(lines[i+2])
						f.write(lines[i+3])
						f.write(lines[i+4])
						f.write(lines[i+5])
						i += 6
						#Writing Regulators information :
						#--------------------------------
						f.write("\n")
						while (("M2       S") not in lines[i] and ("M2       U") not in lines[i]) and i<len(lines)-1:
							if lines[i].startswith("R " + MachineName):
								f.write(lines[i])
								f.write(lines[i+1])
								f.write(" \n \n \n")
							i += 1

	return generators

def main():
	"""
		From a list of eurostag network situation this function extract all the 
		generators that are used in these simulations.  
	"""

	generators = np.zeros(shape=(0,2)) # The list of generators

	# Creating folder where the results are stored :
	# ----------------------------------------------
	if not os.path.exists('rez'):
		os.makedirs('rez')

	# Preparing sim.dta file :
	# ------------------------
	with open('./rez/sim.dta','w') as file_sim:
		file_sim.write("HEADER     "+ time.strftime("%d/%m/%y") +" 5.1\n \n")

	for file_dta in glob.glob('*.dta'):
		
		file_dict = FfindDict(file_dta)
		generators = FreadWritedtaFile(file_dta, file_dict, generators)
		
		time.sleep(0.001)
		print "File is Done !"

	print 'sim.dta has been written'

	np.savetxt('./rez/generators.txt',generators, fmt='%s')
	print 'generators files have been written'

if __name__=='__main__':
	main()
