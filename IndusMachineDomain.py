# -*- coding: utf-8 -*-

# Copyright (c) 2016, Pierre Saikaly  (saikalypierre@gmail.com)
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#===========================#
# created on 06 june 2016
#===========================#

# Import Python dependencies :
# ----------------------------
import imp
import os
import sys
import numpy as np
from numpy import linalg
import math as m
from math import pi
from math import sqrt
from math import fabs
import random as rand
import time
import shutil
from operator import attrgetter
import subprocess
import multiprocessing as mp

def FwriteFiles(MachinesdtaPath,generators):
	""" 
		Create a file for each machine in reference file :
		Inputs  :
			- MachinesdtaPath : location of reference file
		Outputs : 
			- n single machine dta files
		Used in :
			- main
	"""
	regulators = []       # list of regulators for the generator
	generators_in  = np.zeros(shape=(0,4))   # generators that are in the MachinesdtaPath
	generators_out = np.zeros(shape=(0,2))   # generators that are out of MachinesdtaPath

	#Opening Reference file :
	#------------------------
	with open(MachinesdtaPath) as Machinefile:
		lines = Machinefile.readlines()
		#Looking for machines :
		#----------------------
		for i, line in enumerate(lines):
			if line.startswith("M2       U") or line.startswith("M2       S"):
				#Correcting node data :
				#----------------------
				node = lines[i+1].split(" ", 2)
				node[1] = "N1".ljust(len(node[1]))
				lines[i+1] = " ".join(node)
				#Writing specific machine file :
				#--------------------------------
				MachineName = lines[i+1].split(None, 1)[0]
				if MachineName in generators[:,1]:
					gen_in = np.append(generators[MachineName==generators[:,1],:],lines[i+1].split()[2])
					gen_in = np.append(gen_in,lines[i+1].split()[3])
					generators_in = np.vstack((generators_in,gen_in))
					MachineRefName = generators[MachineName==generators[:,1],0][0]
					if not os.path.exists(MachineRefName):
						os.makedirs(MachineRefName)

						time.sleep(0.1)

						filename = 'simTest' + '.dta'
						filename = os.path.join(MachineRefName, filename)
						with open(filename, 'w') as f:
							#Writing Header :
							#----------------
							f.write("HEADER     "+ time.strftime("%d/%m/%y") +" 5.1\n \n")
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
									regulators.append(lines[i+1].split()[0])
									f.write(" \n \n \n")
								i += 1
								time.sleep(0.01)
							#Writing Network information :
							#-----------------------------
							f.write("I1 \nN2             0.     0.02      90.     100. \n \n \nLOADP   1\n         1              1.       1. \n \n \n \n \n"
								"CH \n1        W\n \n \n")

						FcopyRegulators(regulators,MachineRefName)
						FwriteSeq(MachineRefName)
						time.sleep(0.1)
						regulators = []
						print MachineName, " is done !"
					else:
						print MachineName, " already OK"
				else:
					generators_out = np.vstack((generators_out,generators[MachineName==generators[:,1],:]))

	return generators_in, generators_out

def FreadGenerators(GeneratorsPath):
	""" 
		Create a file for each machine in reference file :
		Inputs  :
			- GeneratorsPath : location of generators to process
		Outputs : 
			- Generators : array of generators with dictionary and internal name
		Used in :
			- main
	"""	
	generators = np.zeros(shape=(0,2))

	with open(GeneratorsPath,'r') as GeneratorsFile:
		lines = GeneratorsFile.readlines()

	for line in lines:
		generators = np.vstack((generators,line.split()))

	return generators 

def FcopyRegulators(regulators, MachineName):
	""" 
		Copy the regulators specified :
		Inputs  :
			- regulators : list of regulators to copy
		Outputs : 
			- None, the regulators are copied in the corresponding folder of the machine
		Used in :
			- FwriteFiles
	"""
	for regulator in regulators:
		for basename in os.listdir('.'):
			if basename.startswith(regulator.lower()):
				pathname = os.path.join('.', basename)
				if os.path.isfile(pathname):
					shutil.copy(pathname, MachineName)
		time.sleep(0.01)

def FwriteSeq(MachineName):
	""" 
		Create the seq file for Eurostag simulation
		Inputs  :
			- MachineName
		Outputs : 
			- None, the seq file is created in the folder of the machine
		Used in :
			- 
	"""
	filename = os.path.join(MachineName,'simTest.seq')

	with open(filename,'w') as f:
		f.write("HEADER     "+ time.strftime("%d/%m/%y") +" 5.1\n \n")
		f.write('PARAM\n')
		f.write('         0  0  1\n')
		f.write('                 0.0001         0.0001          0.001             0.             0. 1\n')
		f.write('\n')
		f.write('TIME\n')
		f.write('               0.000001            60.             0.\n')
		f.write('\n')
		f.write('EVENTS\n')
		f.write('   1000. STOP')
		f.write(' \n \n')

def FcallProgram(generator):
	""" 
		Call MachineStableDomainSearch.py for the generator in input 
		Inputs  :
			- generator
		Outputs : 
			- outputs of MachineStableDomainSearch.py
		Used in :
			- main
	"""

	print "Currently working on : ", generator[0]
	# Moving to file of generator :
	# -----------------------------
	os.chdir(generator[0])

	# Lauching python process :
	# -------------------------
	python_process = os.path.join("..","MachineStableDomainSearch.py -v 4 -i 0")
	time_start = time.clock()
	with open("output.txt",'w') as f:
		try:
			subprocess.call("python " + python_process, stdout=f, stderr=subprocess.STDOUT)
		except:
			f.write("WARNING : La recherche n'a pas pu être lancé")
	time_elapsed = (time.clock() - time_start)

	print generator[0], " is Done : ", time_elapsed 
	time.sleep(0.1)
	# Moving back to parent folder
	os.chdir("..")

def FappendResults(generators):
	"""
		Append the results of computed domain to a single text file
		Inputs : 
			- generators : list of inpur generators
		Ouputs : 
			- ampl_generators_domains.txt 
		Used in :
			- main
	"""

	generator_prob = np.zeros(shape=(0,4)) # List of generators with problems

	with open("ampl_generators_domains.txt", "w") as f: 
		f.write("# " + time.strftime("%d/%m/%y") + "\n")
		f.write("#num id P(MW) Q(MVar) V(kV) RHS(lt) " + "Vnominal(kV)" + " " + "id_internal" + "\n")

	with open("ampl_generators_domains.txt", "a") as f: 
		for generator in generators:
			resultfile = os.path.join(generator[0],"ampl_generators_domains.txt")
			try:
				with open(resultfile,'r') as f_res:
					lines = f_res.readlines()
					for line in lines:
						if '#' not in line:
							f.write(line)
			except IOError:
				generator_prob = np.vstack((generator_prob,generator))

			time.sleep(0.1)

	return generator_prob
if __name__ == '__main__':

	# Get list of worker from cpu count
	mp.freeze_support()
	PROCESSES = mp.cpu_count()-1

	RefdtaFile     = sys.argv[1]  
	GeneratorsFile = sys.argv[2]  
	DictFile       = sys.argv[3]  	
	RefdtaPath     = os.path.join(RefdtaFile)     # Path of dtaFile containing machine to test
	DictPath       = os.path.join(DictFile)       # Path of file containing the generators
	GeneratorsPath = os.path.join(GeneratorsFile) # Path of file containing generators to test

	generators    = FreadGenerators(GeneratorsPath)
	generators_in, generators_out = FwriteFiles(RefdtaPath, generators)

	np.savetxt('generators_in_info.txt', generators_in, fmt='%s')

	print "Files have all been written \n"

	pool = mp.Pool(processes=PROCESSES)
	for generator in generators_in:
		pool.apply_async(FcallProgram, args=(generator,))
		time.sleep(0.1)
	pool.close()
	pool.join()

	# waiting for processes to close
	time.sleep(1)
	
	print "Fin du calcul - Concatenation des resultats"
	generators_prob = FappendResults(generators_in)

	print "Liste des generateurs ayant eu un prb : \n" , generators_prob
	np.savetxt('generators_pb.txt',generators_prob, fmt='%s')
	print "Un export a ete fait dans le fichier :\n > generators_pb.txt \n"
	print "Liste des generateurs non inclu dans le .dta : \n" , generators_out
	np.savetxt('generators_out.txt',generators_out, fmt='%s')
	print "Un export a ete fait dans le fichier :\n > generators_out.txt \n"


	