
# -*- coding: utf-8 -*-

# Copyright (c) 2016, Pierre Saikaly  (saikalypierre@gmail.com)

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this file, 
# You can obtain one at http://mozilla.org/MPL/2.0/.

#===========================#
# created on 12 july 2016
#===========================#

#Import Python dependencies :
#----------------------------
import os
import sys
import numpy as np
import random as rand
from numpy import linalg
import math as m
from math import pi
from math import sqrt
from math import fabs
import time
from scipy.spatial import ConvexHull
import argparse

# Import custom dependencies :
# ----------------------------

# Import Python visualisation dependencies :
# -----------------------------------------
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def FreadMaxMin(dataPath):
	"""
		Open computed domain file and exctract max/min active and reactive power 
		Input :
			- dataPath : path of domain file
		Outputs : 
			- min/max P
			- min/max Q
		Used in :
			- 
	"""

	# Reading data file for min and max :
	# -----------------------------------
	with open(dataPath, 'r') as datafile:
		lines = datafile.readlines()

	# /!\ if ampl_gen.. structure is changed the following is no longer true /!\
	lines[1] = lines[1].translate(None, '[]')
	lines[2] = lines[2].translate(None, '[]')
	lines[3] = lines[3].translate(None, '[]')
	lines[4] = lines[4].translate(None, '[]')
	lines[5] = lines[5].translate(None, '[]')
	lines[6] = lines[6].translate(None, '[]')

	Pmax = float(lines[1].split()[3])
	Pmin = float(lines[2].split()[3]) 
	Qmax = float(lines[3].split()[4])
	Qmin = float(lines[4].split()[4])
	Umax = float(lines[5].split()[5])
	Umin = float(lines[6].split()[5])

	# Enlarging field for better view :
	# ---------------------------------
	Qmax *= 1.2
	Qmin = Qmin - 0.2*abs(Qmin) 
	Umax *= 1.2
	Umin *= 0.8

	return Pmax, Pmin, Qmax, Qmin, Umax, Umin

def FreadDiagramUQ(diagUQPath,typediag):
	"""
		Extract domain information from reference excel file and store the points extracted
		If the reference file has a different structure this function will not work
		Inputs :
			- diagUQPath : path of the reference file
			- typediag   : Extract data for the specified type
		Outputs : 
			- data   : points written in file
			- P_draw : active power for which the points are computed
		Used in :
			- 
	"""
	convert = False     # indicator if data needs to be converted
	aux     = 1         # indicator if aux conversion is possible

	with open(diagUQPath, 'r') as diagUQfile:
		lines = diagUQfile.readlines()

	# Exctracting data :
	# ------------------
	for i, line in enumerate(lines):
		if line.startswith('VALPTRACE'):
			P_draw = float(line.split(';')[2])
		if (line.startswith('LTRACE') and line.split(';')[2]=='R\n'):
			convert = True

		if (convert == True):
			if line.startswith('PAUX'):
				try:
					Paux = float(line.split(';')[2].strip())
				except ValueError :
					Paux = 0
					aux = 0
			if line.startswith('QAUX'):
				try:
					Qaux = float(line.split(';')[2].strip())
				except ValueError :
					Qaux = 0
					aux = 0
			if line.startswith('XTFO'):
				try:
					XTFO = float(line.split(';')[2].strip())
				except ValueError :
					XTFO = 0
			if line.startswith('UBASE'):
				try:
					Ubase = float(line.split(';')[2].strip())
				except ValueError :
					Ubase = 0
			if line.startswith('SBASE'):
				try:
					Sbase = float(line.split(';')[2].strip())
				except ValueError :
					Sbase = 0
			if line.startswith('XLIG'):
				try:
					Xlig = float(line.split(';')[2].strip())
				except ValueError :
					Xlig = 0
			if line.startswith('RLIG'):
				try:
					Rlig = float(line.split(';')[2].strip())
				except ValueError :
					Rlig = 0
			if line.startswith('BLIG'):
				try:
					Blig = float(line.split(';')[2].strip())*10**(-6)
				except ValueError :
					Blig = 0
			if line.startswith('USN'):
				try:
					Usn = float(line.split(';')[2].strip())
				except ValueError :
					Usn = 0

		if line.startswith('TYPEDIAG'):
			typesdiag = line.split(';')
			for j in range(2,len(typesdiag)):
				if typediag == typesdiag[j].strip():
					position = j-2
					break

		if line.startswith('NBPTS'):
			if not line.split(';')[2 + position]:
				data_size = 0
			else:
				data_size = int(line.split(';')[2 + position])
				data = np.zeros(shape=(data_size,2))
				
		if line.startswith('UNITE'):
			if (data_size != 0):
				for k in range(0,data_size):
					data[k,:] = [float(lines[i+k+1].split(';')[2 + position]), float(lines[i+k+1].split(';')[3 + position])]
				break
			else:
				data = np.zeros(shape=(0,2))
				k = 0
				while 'FIN' not in lines[i+k+1]:
					data = np.vstack((data,[float(lines[i+k+1].split(';')[2 + position]), float(lines[i+k+1].split(';')[3 + position])]))
					k += 1

	# If data is given at the network side it needs to be converted to power generator side
	if (convert == True):
		data = FcalcNetworktoGen(data, P_draw, Paux, Qaux, XTFO, Ubase, Sbase, Xlig, Rlig, Blig, Usn, aux)

	return data, P_draw

def FcalcNetworktoGen(data, P_draw, Paux, Qaux, XTFO, Ubase, Sbase, Xlig, Rlig, Blig, Usn, aux):
	"""
		Convert data from the network side to the generator side.
		This is necessary because the computed domain are computed on the generator side
		Inputs : 
			- data   : points in the reference diagram
			- P_draw : active power 
			- caracteristics of the generators
		Outputs : 
			- data_cv : converted data points 
		Used in :
			- FreadDiagramUQ
	"""

	data_int = np.zeros(shape=(len(data),2))   # calculation variable
	data_cv  = np.zeros(shape=(len(data),2))   # output
	Pint     =  np.zeros(shape=(len(data),1))  # calculation variable
	Pcv      = np.zeros(shape=(len(data),1))   # calculation variable
	
	for i in range(0,len(data)):
		Intensite1 = sqrt(P_draw**2 + data[i,1]**2) / data[i,0]
		Pint[i] = P_draw + Rlig*Intensite1**2
		data_int[i,1] = data[i,1] + Xlig*Intensite1**2 - 2*Blig*data[i,0]**2
		data_int[i,0] = sqrt(Pint[i]**2 + data_int[i,1]**2)/Intensite1
	
	if (aux == 1):
		for i in range(0,len(data)):
			Intensite2   = sqrt((Pint[i] + Paux)**2 + (data_int[i,1] + Qaux)**2)/data_int[i,0]*Ubase/Usn
			Pcv[i] = Pint[i] + Paux
			data_cv[i,1] = data_int[i,1] + Qaux + Intensite2**2*XTFO/100*Usn**2/Sbase
			data_cv[i,0] = sqrt(Pcv[i]**2 + data_cv[i,1]**2)/Intensite2
	else:
		for i in range(0,len(data)):
			Intensite2   = sqrt((Pint[i])**2 + (data_int[i,1])**2)/data_int[i,0]*Ubase/Usn
			Pcv[i] = Pint[i] + Paux
			data_cv[i,1] = data_int[i,1] + Qaux + Intensite2**2*XTFO/100*Usn**2/Sbase
			data_cv[i,0] = sqrt(Pint[i]**2 + (data_cv[i,1]-Qaux)**2)/Intensite2

	return data_cv

def FinDomain(X,data_raw,P_draw):
	"""
		Check if the point X is inside the domain made with data_raw at P_draw
		Inputs : 
			- X : point to be tested
			- data_raw : constraints 
			- P_draw : Active power at which the test is made
		Outputs :
			- True/False
		Used in :

	"""	
	# Local Variable :
	# ----------------
	eps = 0.001    # precision

	# Testing point : 
	# ---------------
	for i in range(0,len(data_raw)):
		if (np.dot(data_raw[i,1:3],X) - data_raw[i,3] + data_raw[i,0]*P_draw >= eps):
			return False
			break
	
	return True

def FcomputeUQdiag(data_raw, P_draw):
	"""
		Compute diagram UQ from the constrains at the specified active power.
		it return the points sorted in order to visualise them
		Inputs : 
			- data_raw : constrains
			- P_draw   : Active power 
		Outputs : 
			- Points inside the domain sorted 
		Used in :
			- 
	"""
	eps    = 0.00000000001   # tolerance for system
	eps2   = 0.000001        # tolerance for duplicate
	points = np.zeros(shape=(0,2))

	# Computing intersections points :
	# --------------------------------
	for i in range(0,len(data_raw)):
		for j in range(0,len(data_raw)):
			M_sys = np.vstack((data_raw[i,1:3],data_raw[j,1:3]))
			delta = np.linalg.det(M_sys)
			if (abs(delta) > eps):
				X = np.linalg.solve(M_sys,[data_raw[i,3] - data_raw[i,0]*P_draw, data_raw[j,3] - data_raw[j,0]*P_draw])
				# Check if X is in domain : 
				# -------------------------
				if (FinDomain(X, data_raw, P_draw)):
					doublon = False
					for l in range(0,len(points)):
						if np.linalg.norm(X-points[l,:])<eps2:
							doublon = True
							break
					if(doublon!=True):
						points = np.vstack((points,X))

	# Computing Hull of points to order them :
	# ----------------------------------------
	Hull = ConvexHull(points)
	points_sorted = points[Hull.vertices,:]

	return points_sorted

def main(diagUQPath, typediag, drawPlane, domainsupport, domaincorners, domainrandom, domainref):
	"""
		Compute and display the UQ Diagram with several options and types of vizualisation
		Inputs :
			- diagUQPath : The reference UQ diagram comes in a .csv folder. 
			It needs to have a specific syntax to be used.
			- typeidag   : the type of data that is going to be extracted from 
			the reference UQ file the options are :
			"ZEC RPT", "ZFN RPT RTE", "ZEC RST", "ZFN RST RTE", "ZFN"
			- domain type : different options for visualisations :
			0 : Support points UQ diagram
			1 : Corner points UQ diagram
			2 : Support points + Corner points UQ diagram
			3 : Support points + Corners points + Sample points UQ diagram
			4 : Ref + Support points + Corners points + Sample points UQ diagram
			- drawPlane : draw the projection of all the planes in the UQ diagram 
			Use only when there is one diagram else it becomes a mess
		Outputs : 
			- figures 
	"""

	# Setting up Visualisation :
	# ---------------------------
	fig, ax = plt.subplots(figsize=(9, 9))
	fig.subplots_adjust(left = 0.1, bottom = 0.1,
							right = 0.9, top = 0.9, wspace = 0, hspace = 0.01)

	dataPath = "ampl_generators_domains.txt"

	# Reading reference UQ diagram and extracting power to draw :
	# -----------------------------------------------------------
	dataUQRef, P_draw = FreadDiagramUQ(diagUQPath,typediag)

	Pmax, Pmin, Qmax, Qmin, Umax, Umin = FreadMaxMin(dataPath)

	if (domainref == 1):
		P_draw_array = [P_draw]
	else:
		P_draw_array = np.linspace(Pmin+0.1*abs(Pmin), 0.9*Pmax, num = 15)

	for i, P_draw in enumerate(P_draw_array):
		if (len(P_draw_array) > 1):
			plt.subplot(5,3,i+1)

		print P_draw
						
		if (domaincorners == 1):
			# Draw UQ diagram from corner points
			dataPath = "ampl_generators_domains_coins.txt"

			# Reading data
			data_raw = np.loadtxt(dataPath, comments='#', usecols=(2,3,4,5,6))

			# Computing UQ diagram from data for the selected power :
			points = FcomputeUQdiag(data_raw, P_draw)

			# Creating visualisation :
			# ------------------------
			U = np.linspace(Umin, Umax)
			plt.fill(points[:,1], points[:,0], color='teal', alpha=1, label='Tangent')
			# Draws all projected constraints
			if (drawPlane == 1):
				for i in range(0,len(data_raw)):
					if (abs(data_raw[i,1]) > 0.00000001):
						plt.plot(U, (data_raw[i,3] - data_raw[i,0]*P_draw - data_raw[i,2]*U)/data_raw[i,1], lw=0.5, color='black')

		if (domainsupport == 1):
			# Draw UQ diagram from support points
			dataPath = "ampl_generators_domains.txt"

			# Reading data
			data_raw = np.loadtxt(dataPath, comments='#', usecols=(2,3,4,5,6))

			# Computing UQ diagram from data for the selected power :
			points = FcomputeUQdiag(data_raw, P_draw)

			# Creating visualisation :
			# ------------------------
			U = np.linspace(Umin, Umax)
			plt.fill(points[:,1], points[:,0], color='purple', alpha=1, label='Support')
			# Draws all projected constraints
			if (drawPlane == 1):
				for i in range(0,len(data_raw)):
					if (abs(data_raw[i,1]) > 0.00000001):
						plt.plot(U, (data_raw[i,3] - data_raw[i,0]*P_draw - data_raw[i,2]*U)/data_raw[i,1], lw=0.5, color='black')

		if (domainrandom == 1):
			# Draw UQ Diagram from sample points
			dataPath = "ampl_generators_domains_MC.txt"

			# Reading data
			data_raw = np.loadtxt(dataPath, comments='#', usecols=(2,3,4,5,6))

			# Computing UQ diagram from data for the selected power :
			points = FcomputeUQdiag(data_raw, P_draw)

			# Creating visualisation :
			# ------------------------
			U = np.linspace(Umin, Umax)
			plt.fill(points[:,1], points[:,0], color='yellow', alpha=1, label='Random ')
			# Draws all projected constraints
			if (drawPlane == 1):
				for i in range(0,len(data_raw)):
					if (abs(data_raw[i,1]) > 0.00000001):
						plt.plot(U, (data_raw[i,3] - data_raw[i,0]*P_draw - data_raw[i,2]*U)/data_raw[i,1], lw=0.5, color='black')

		# Drawing reference UQ diagram :
		# ------------------------------
		if (domainref == 1):
			# for i in range(0,len(dataUQRef)):
			# 	plt.scatter(dataUQRef[i,0],dataUQRef[i,1], color='green', s=10)
			plt.fill(dataUQRef[:,0],dataUQRef[:,1], color='green', alpha=1, label='Reference')
			plot_title = "Diagramme UQ a " + str(round(P_draw,1)) + " MW"
			plt.xlabel('U Tension (kV)', fontsize='x-large')
			plt.ylabel('Q Puissance reactive (MVar)', fontsize='x-large')

		# Wrapping up visualisation :
		# ---------------------------
		plt.axis([Umin, Umax, Qmin, Qmax])
		if (domainref == 0):
			plt.xlabel('U', fontsize='small')
			plt.ylabel('Q ', fontsize='small')
			# plt.title(plot_title, fontsize='xx-large')
			# ax.axis('off')
			# plt.axis('off')
			
		plt.draw()
		# plt.show()
		
	plt.legend(loc='best', shadow=True, fontsize='large')
	if (len(P_draw_array) > 1):
		fig.savefig("Diag_UQ_multiple.png")
	else:
		plt.title(plot_title, fontsize='xx-large')
		fig.savefig("Diag_ref.png")


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Create a plane projection, you must be first be inside the simulation directory and run "python ..\CreatePlaneDomain.py + arguments"')
	parser.add_argument('-p','--UQpath', help='Path of the reference diagram',
			 required=True)
	parser.add_argument('-t','--type', help='Type of diagram (* the one commonly used) : either "ZEC RPT", *"ZFN RPT RTE", "ZEC RST", "ZFN RST RTE", "ZFE"',
			 required=True)
	parser.add_argument('-d','--dplane', help='Draw the projection of all the planes : 0 off, 1 on',required=True)
	parser.add_argument('-s','--support', help='Draw the supports projection : 0 off, 1 on',required=True)
	parser.add_argument('-c','--corners', help='Draw the corners projection : 0 off, 1 on',required=True)
	parser.add_argument('-a','--random', help='Draw the random projection : 0 off, 1 on',required=True)
	parser.add_argument('-r','--reference', help='Draw the reference diagram : 0 off, 1 on, if on only plot for the power of diagram',required=True)
	args = parser.parse_args()


	# Reading input arguments :
	# -------------------------
	diagUQPath = args.UQpath
	typediag   = args.type
	drawPlane  = int(args.dplane)
	domainsupport = int(args.support)
	domaincorners = int(args.corners)
	domainrandom  = int(args.random)
	domainref     = int(args.reference)

	main(diagUQPath, typediag, drawPlane, domainsupport, domaincorners, domainrandom, domainref)






















