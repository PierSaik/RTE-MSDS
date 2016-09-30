# -*- coding: utf-8 -*-

# Copyright (c) 2016, Pierre Saikaly  (saikalypierre@gmail.com)
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#===========================#
# created on 14 april 2016
#===========================#

# Import Python dependencies :
# ----------------------------
import imp
import os
import sys
import shutil
import numpy as np
from numpy import linalg
import math as m
from math import pi
from math import sqrt
from math import fabs
import random as rand
import time
from operator import attrgetter
import scipy.optimize
import re
import subprocess

# Import Python visualisation libraries :
# ---------------------------------------
# import matplotlib as mpl
# mpl.use('Agg')
# import matplotlib.pyplot as plt
# import matplotlib.cm as cm
# from mpl_toolkits.mplot3d import Axes3D
# from itertools import product, combinations
# from matplotlib.patches import FancyArrowPatch
# from mpl_toolkits.mplot3d import proj3d
# from matplotlib.path import Path


from scipy.spatial import ConvexHull


# Import custom dependencies :
# -----------------------------
import CreateMachineDomain

# Defining vector visualisation class :
# -------------------------------------
# class Arrow3D(FancyArrowPatch):
#     def __init__(self, xs, ys, zs, *args, **kwargs):
#         FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
#         self._verts3d = xs, ys, zs

#     def draw(self, renderer):
#         xs3d, ys3d, zs3d = self._verts3d
#         xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
#         self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
#         FancyArrowPatch.draw(self, renderer)

# Data class :
# ------------
class PointBoundary():
	"""Define point of boundary in pu and si with its normal"""

	def __init__(self,p_input,n_input,unit,P_boundaries,Q_boundaries,U_boundaries):
		"""
			Define point and normal in both pu and si :
				Inputs :
					p_input : P,Q,U point
					n_input : P,Q,U vector
					unit    : si or pu
				Output :
					point_si : P,Q,U in si
					point_pu : P,Q,U in pu
					normal_si : P,Q,U in si
					normal_si : P,Q,U in pu
		"""
		if unit=="pu":
			self.point_pu = p_input
			self.normal_pu = n_input
			self.point_si = Conversion_pu_si(self.point_pu,P_boundaries,Q_boundaries,U_boundaries)
			self.normal_si = Conversion_pu_si(self.normal_pu,P_boundaries,Q_boundaries,U_boundaries)
		elif unit=="si":
			self.point_si = p_input
			self.normal_si = n_input
			self.point_pu = Conversion_si_pu(self.point_si,P_boundaries,Q_boundaries,U_boundaries)
			self.normal_pu = Conversion_si_pu(self.normal_si,P_boundaries,Q_boundaries,U_boundaries)
		else:
			print "Specifier l'unite : pu ou si"

class PointDomain():
	"""Define point of domain in pu and si"""

	def __init__(self,p_input,unit,P_boundaries,Q_boundaries,U_boundaries):
		"""
			Define point in both pu and si :
				Inputs :
					p_input : P,Q,U point
					unit    : si or pu
				Output :
					point_si : P,Q,U in si
					point_pu : P,Q,U in pu
		"""
		if unit=="pu":
			self.point_pu = p_input
			self.point_si = Conversion_pu_si(self.point_pu,P_boundaries,Q_boundaries,U_boundaries)
		elif unit=="si":
			self.point_si = p_input
			self.point_pu = Conversion_si_pu(self.point_si,P_boundaries,Q_boundaries,U_boundaries)
		else:
			print "Specifier l'unite : pu ou si"

def Conversion_pu_si(data_pu,P_boundaries,Q_boundaries,U_boundaries):
	"""
		Converts data from pu to si :
		Inputs :
			- data_pu : P,Q,U data in pu
			- U_nominal : nominal tension
			- P_nominal : nominal power
		Output :
			- data_si : Converted point in si
		Used in :
			- 
	"""
	# Creating output :
	# -----------------
	data_si = np.zeros(3)
	# Converting :
	# ------------
	data_si[0] = (data_pu[0]) * (P_boundaries[1] - P_boundaries[0]) + P_boundaries[0]
	data_si[1] = (data_pu[1]) * (Q_boundaries[1] - Q_boundaries[0]) + Q_boundaries[0]
	data_si[2] = (data_pu[2]) * (U_boundaries[1] - U_boundaries[0]) + U_boundaries[0]

	return data_si

def Conversion_si_pu(data_si,P_boundaries,Q_boundaries,U_boundaries):
	"""
		Converts data from pu to si :
		Inputs :
			- data_si : P,Q,U data in si
			- U_nominal : nominal tension
			- P_nominal : nominal power
		Output :
			- data_pu : Converted point in pu
		Used in :
			- 
	"""	
	# Creating output :
	# -----------------
	data_pu = np.zeros(3)
	# Converting :
	# ------------
	data_pu[0] = ((data_si[0] - P_boundaries[0]) / (P_boundaries[1] - P_boundaries[0])) 
	data_pu[1] = ((data_si[1] - Q_boundaries[0]) / (Q_boundaries[1] - Q_boundaries[0])) 
	data_pu[2] = ((data_si[2] - U_boundaries[0]) / (U_boundaries[1] - U_boundaries[0])) 

	return data_pu

# Plane fitting from n points :
# -----------------------------
def fitPLaneLTSQ(XYZ):
    # Fits a plane to a point cloud, 
    # Where Z = aX + bY + c        ----Eqn #1
    # Rearanging Eqn1: aX + bY -Z +c =0
    # Gives normal (a,b,-1)
    # Normal = (a,b,-1)
    [rows,cols] = XYZ.shape
    G = np.ones((rows,3))
    G[:,0] = XYZ[:,0]  #X
    G[:,1] = XYZ[:,1]  #Y
    Z = XYZ[:,2]
    (a,b,c),resid,rank,s = np.linalg.lstsq(G,Z) 
    normal = (a,b,-1)
    nn = np.linalg.norm(normal)
    normal = normal / nn
    return normal

def fitPlaneOptimize(XYZ):
    def residiuals(parameter,f,x,y):
        return [(f[i] - model(parameter,x[i],y[i])) for i in range(len(f))]

    def model(parameter, x, y):
        a, b, c = parameter
        return a*x + b*y + c

    X = XYZ[:,0]
    Y = XYZ[:,1]
    Z = XYZ[:,2]
    p0 = [1., 1.,1.] # initial guess
    result = scipy.optimize.leastsq(residiuals, p0, args=(Z,X,Y))[0]
    normal = result[0:3]
    nn = np.linalg.norm(normal)
    normal = normal / nn
    return normal

def fitPlaneSVD(XYZ):
    [rows,cols] = XYZ.shape
    # Set up constraint equations of the form  AB = 0,
    # where B is a column vector of the plane coefficients
    # in the form b(1)*X + b(2)*Y +b(3)*Z + b(4) = 0.
    p = (np.ones((rows,1)))
    AB = np.hstack([XYZ,p])
    [u, d, v] = np.linalg.svd(AB,0)        
    B = v[3,:];                    # Solution is last column of v.
    nn = np.linalg.norm(B[0:3])
    B = B / nn
    return B[0:3]

def fitPlaneSolve(XYZ):
    X = XYZ[:,0]
    Y = XYZ[:,1]
    Z = XYZ[:,2] 
    npts = len(X)
    A = np.array([ [sum(X*X), sum(X*Y), sum(X)],
                   [sum(X*Y), sum(Y*Y), sum(Y)],
                   [sum(X),   sum(Y), npts] ])
    B = np.array([ [sum(X*Z), sum(Y*Z), sum(Z)] ])
    normal = np.linalg.solve(A,B.T)
    nn = np.linalg.norm(normal)
    normal = normal / nn
    return normal.ravel()

def fitPlaneEigen(XYZ):
	average=sum(XYZ)/XYZ.shape[0]
	covariant=np.cov(XYZ - average)
	eigenvalues,eigenvectors = np.linalg.eig(covariant)
	want_max = eigenvectors[:,eigenvalues.argmax()]
	(c,a,b) = want_max[3:6]
	normal = np.array([a,b,c])
	nn = np.linalg.norm(normal)
	return normal / nn  

# Program functions :
# -------------------	
def FdynamoSimu(parametersPath,X, U_nominal):
	""" Call Dynamo simulator
		Inputs :
			- X       : machine parameters [P,Q,U]
			- seqPath : Machine .seq file 
			- dtaPath : Machine .dta file
			- echPath : Machine .ech file
		Output : 
			- alpha : Eurostag status
			- X1    : Eurostag stable point

		Used in :
			- Fexploration
			- Fraffinement
	"""
	# # Removing previous
	# try:
	# 	shutil.rmtree('outputs')
	# except:
	# 	pass

	dynamo_output = np.zeros(shape=(1,3))
	outputsPath = os.path.join(".","outputs","debug","globalInit","dumpInitValues-MACHINE.txt")

	FeditParameters(parametersPath, X, U_nominal)

	# call dynamo :
	# -------------
	with open("dynamo_output.txt", 'w') as f:
		subprocess.call("export DYNAMO_INSTALL_DIR=~/dynamo/", shell=True, stdout=f, stderr=subprocess.STDOUT)
	with open("dynamo_output.txt", 'a') as f:
		subprocess.call(["./bin/execDynamo.sh", "jobs", "noIIDM.jobs"], stdout=f, stderr=subprocess.STDOUT)

	# Reading outputs
	# status, var_instable = FreadOutputs(outputsPath)

	# print "Max variable instalbe : ", max(var_instable[:,1])

	return X

def FeditParameters(parametersPath,X,U_nominal):
	""" 
		Edit parameters.par file :
		Inputs :
			- parametersPath : path to parameters.par file
			- X : [P,Q,U] point in [MW,MVar,kV]
		Output : 
			- parameters.par file with [P,Q,U,theta] point in per unit

		Used in :
			- 
	"""
	# Internal variable :
	# -------------------
	Snref = 100. # From MW, Mvar to p.u

	# Parameters name :
	# -----------------
	name_P  = 'M2S_P0'
	name_Pc = 'TRANSFO_P10'
	name_Q  = 'M2S_Q0'
	name_Qc = 'TRANSFO_Q10'
	name_U  = 'M2S_U0'
	name_Uc = 'TRANSFO_U10'
	name_U2 = 'UVGT_U0'

	# Regex :
	# -------
	reg_name  = re.compile(r'(?<=name=")([^">]+)')
	reg_value = re.compile(r'(?<=value=")([^">]+)')

	with open(parametersPath,'r') as parametersFile_in:
		lines = parametersFile_in.readlines()

	for i, line in enumerate(lines):
		match = reg_name.search(str(line))
		if match != None:
			if ((match.group() == name_P) or (match.group() == name_Pc)):
				lines[i] = reg_value.sub(str(X[0]/Snref),str(line))
			if ((match.group() == name_Q) or (match.group() == name_Qc)):
				lines[i] = reg_value.sub(str(X[1]/Snref),str(line))
			if ((match.group() == name_U) or (match.group() == name_Uc) or (match.group() == name_U2)):
				lines[i] = reg_value.sub(str(X[2]/U_nominal),str(line))

	with open(parametersPath,'w') as parametersFile_out:
		for line in lines:
			parametersFile_out.write(line)

def FreadOutputs(outputsPath):
	""" 
		Read outputs file :
		Inputs :
			- outputsPath : 
		Output : 
			- 
		Used in :
			- 
	"""
	# Parameters :
	# ------------
	epsilon = 0.00009
	status  = 0
	var_instable = np.zeros(shape=(0,2))

	with open(outputsPath,'r') as outputsFile:
		lines = outputsFile.readlines()

	for i, line in enumerate(lines):
		if ('====== INIT VARIABLES VALUES ======' in line):
			j = 1
			while ('====== INIT DISCRETE VARIABLES VALUES ======' not in lines[i+j]):
				varia = lines[i+j].split()[0] 
				value = float(lines[i+j].split()[7]) 
				if (abs(value) > epsilon):
					status = 1
					var_instable = np.vstack((var_instable,[varia,value])) 
				j += 1

	return status, var_instable

def Fexploration(O, B, Esp, parametersPath, U_nominal, P_nominal, P_boundaries, Q_boundaries, U_boundaries):
	""" 
		Compute last stable point in a direction :
		Inputs : 
			- O   : Stable origin
			- B   : Boundary point in direction
			- Esp : Boundary precision
		Output :
			- A   : Last Stable point 

		Used in :
			- FplanTangent
	"""
	# Initialisation des variables internes :
	# ---------------------------------------
	C = np.zeros(3) # point in unstable part
	A = np.zeros(3) # point in stable part
	T = np.zeros(3) # test point
	X = np.zeros(3) # eurostag point result
	d = 0 # distance between stable and unstable point
	outputsPath = os.path.join(".","outputs","debug","globalInit","dumpInitValues-MACHINE.txt")
	# Initialisation :
	# ----------------
	A = O
	C = B
	T = C
	A_pu = Conversion_si_pu(A,P_boundaries,Q_boundaries,U_boundaries)
	C_pu = Conversion_si_pu(C,P_boundaries,Q_boundaries,U_boundaries)
	d = np.linalg.norm(A_pu-C_pu)

	# Last stable point at Eps biais :
	# --------------------------------
	while d>Esp:
		print "Distance point stable/instable : ", d, "\n"
		print "Point test :", T
		# try:
		if True:
			X = FdynamoSimu(parametersPath,T, U_nominal)
			# print "Simulation is OK"
			# time.sleep(1)
		# except:
		# 	print "Simulation is KO"
		# 	time.sleep(1)
		# 	pass

		alpha, var_instable = FreadOutputs(outputsPath)
		shutil.rmtree('outputs')
		print var_instable
		print "alpha =", alpha

		if alpha==0:
			#point is on stable side
			A = T
			T = (C+A)/2
			print ".............. Cote Stable .............\n"
		else:
			#point is on unstable divergent side
			C = T
			T = (C+A)/2
			print ".............. Cote instable .............\n"
		# time.sleep(4)
		A_pu = Conversion_si_pu(A,P_boundaries,Q_boundaries,U_boundaries)
		C_pu = Conversion_si_pu(C,P_boundaries,Q_boundaries,U_boundaries)
		D_pu = C_pu - A_pu
		d = np.linalg.norm(D_pu)

	print "Dernier point stable   : ", A, '\n'
	print "Dernier point instable : ", C, '\n'
	return A

def Fvoisinage(B,P_boundaries,Q_boundaries, U_boundaries, pas, sign1, sign2, U_nominal, P_nominal, precision_vect):
	""" 
		Return points in the vicinity of B
		Inputs :
			- B : point from which the vinicity is computed
			- Boundaries of domain
			- pas : step between points
		Output : 
			- voisins : Point in vicinity of B

		Used in :
			- FplanTangent
	"""	
	voisins    = np.zeros(shape=(3))
	voisins[:] = B

	[pas_P, pas_Q, pas_U] = Conversion_pu_si(precision_vect,P_boundaries,Q_boundaries,U_boundaries) - [P_boundaries[0], Q_boundaries[0], U_boundaries[0]]
	print [pas_P, pas_Q, pas_U]
	if np.any(np.isclose(B[0],P_boundaries)):
		# print "CAS P"
		# time.sleep(1)
		voisins[1] += (1+sign1)*pas_Q*100.
		voisins[2] += (1+sign2)*pas_U*100.

	elif np.any(np.isclose(B[1],Q_boundaries)):
		# print "CAS Q"
		# time.sleep(1)
		voisins[0] += (1+sign1)*pas_P*100.
		voisins[2] += (1+sign2)*pas_U*100.

	elif np.any(np.isclose(B[2],U_boundaries)):
		# print "CAS U"
		# time.sleep(1)
		voisins[0] += (1+sign1)*pas_P*10.
		voisins[1] += (1+sign2)*pas_Q*10.
	else:
		return "Le point ", B, " n'est pas sur le bord du domaine"

	# time.sleep(1)

	return voisins

def FplanTangent(O, theta, phi, pt, P_boundaries, Q_boundaries, U_boundaries, Eps, parametersPath, U_nominal, P_nominal, precision_vect):
	""" 
		Compute tangent plane at specified point/direction
		Inputs :
			- O : stable origin
			- D : Domain min & max boundaries
			- p : number of point allowed for tangent evaluation 
		Output : 
			- A : studied point
			- n : normal vector

		Used in :
			- main
	"""
	#Creating local variables :
	#--------------------------
	A = np.zeros(shape=(pt,3)) #Storage matrix of stable point
	B = np.zeros(shape=(pt,3)) #Storage matrix of unstable point
	results = np.zeros(shape=(0,3))
	M = np.zeros(shape=(3,3)) #Covariance matrix of A
	G = np.zeros(3) #Barycentre
	ecart = Eps*100 # 1/sqrt(2*P_nominal**2+U_nominal**2) # Step to compute other points in vicinity

	method = 2

	# espt  = 0.01*raff #wideness in theta (2*espt)
	# espp  = 0.01*raff #wideness in phi (2*espp)


	# Computing first point on stable domaine boundary :
	# --------------------------------------------------
	B[0,:] = Fdirection(theta,phi,P_boundaries, Q_boundaries, U_boundaries, O, U_nominal, P_nominal)
	A[0,:] = Fexploration(O, B[0,:], Eps, parametersPath, U_nominal, P_nominal, P_boundaries, Q_boundaries, U_boundaries)
	
	results = np.vstack((A[0,:],results))

	if method == 1 :
		# Computing angle step for other points :
		# ---------------------------------------
		distance = np.linalg.norm(Conversion_si_pu(O - B[0,:],P_boundaries,Q_boundaries,U_boundaries))
		epsilon  = m.atan(20 * max(precision_vect[0],precision_vect[1]) / distance)
		epsilonp = m.atan((20 * precision_vect[2] / distance))
		# Computing other directions and stable boundary points :
		# -------------------------------------------------------
		for i in range(1,9):
			signs = "{0:b}".format(i-1).zfill(3) # binary of i with 3 digits
			if int(signs[0])==0:
				sign1 = int(signs[2])*-2
				sign2 = int(signs[1])*-2
			else:
				sign1 = - (int(signs[2]) + int(signs[1]))
				sign2 = int(signs[2]) - int(signs[1]) - 1

			it = 0

			A_pu0 = Conversion_si_pu(A[0,:], P_boundaries, Q_boundaries, U_boundaries)
			A_puT = Conversion_si_pu(A[i,:], P_boundaries, Q_boundaries, U_boundaries)
			dist_voisin = np.linalg.norm(A_pu0 - A_puT)

			while ((dist_voisin>1.1*ecart) or it<1) and it<=5:
				it += 1
				B[i,:]  = Fdirection(theta + (sign1 + 1)*epsilon,phi + (sign2 + 1)*epsilonp,P_boundaries, Q_boundaries, U_boundaries, O, U_nominal, P_nominal)
				A[i,:]  = Fexploration(O, B[i,:], Eps, parametersPath, U_nominal, P_nominal, P_boundaries, Q_boundaries, U_boundaries)
				if it > 1:
					A_pu0 = Conversion_si_pu(A[0,:], P_boundaries, Q_boundaries, U_boundaries)
					A_puT = Conversion_si_pu(A[i,:], P_boundaries, Q_boundaries, U_boundaries)
					dist_voisin = np.linalg.norm(A_pu0 - A_puT)


					if 0.5*ecart>dist_voisin:
						epsilon += epsilon/2
						epsilonp += epsilonp/2
						print "Plus de epsilon"
					elif dist_voisin>1.1*ecart:
						epsilon -= epsilon/2
						epsilonp -= epsilonp/2
						print "Moins de epsilon"

					print "Angles : ", epsilon, epsilonp
					print "Distance entre les points : ", dist_voisin
					print "Critère : ", ecart
					print "pour le point : ", i
					# time.sleep(5)

			if dist_voisin<5*ecart:
				results = np.vstack((results,A[i,:]))

	if method == 2:
		precision_vectin = np.zeros(shape=(3))
		Bbis     = np.zeros(shape=(1,3))
		critere  = Eps * 100
		ecartbis = critere * 5
		Obis = O + 0.70*(O-A[0,:])
		
		for i in range(1,9):
			precision_vectin[:] = precision_vect
			Bbis  = np.zeros(shape=(1,3))
			signs = "{0:b}".format(i-1).zfill(3) # binary of i with 3 digits
			if int(signs[0])==0:
				sign1    = int(signs[2])*-2
				sign2    = int(signs[1])*-2
			else:
				sign1    = - (int(signs[2]) + int(signs[1]))
				sign2    = int(signs[2]) - int(signs[1]) - 1

			it = 0

			A_pu0 = Conversion_si_pu(A[0,:], P_boundaries, Q_boundaries, U_boundaries)
			A_puT = Conversion_si_pu(A[i,:], P_boundaries, Q_boundaries, U_boundaries)
			dist_voisin = np.linalg.norm(A_pu0 - A_puT)
			
			while (( (0.5*critere>dist_voisin or dist_voisin>1.5*critere) and it<=10) or it<1):
				it += 1
				Bbis = Fvoisinage(B[0,:],P_boundaries,Q_boundaries,U_boundaries,ecartbis,sign1,sign2,U_nominal,P_nominal,precision_vectin)
				A[i,:] = Fexploration(O, Bbis, Eps, parametersPath, U_nominal, P_nominal, P_boundaries, Q_boundaries, U_boundaries)

				A_pu0 = Conversion_si_pu(A[0,:], P_boundaries, Q_boundaries, U_boundaries)
				A_puT = Conversion_si_pu(A[i,:], P_boundaries, Q_boundaries, U_boundaries)
				dist_voisin = np.linalg.norm(A_pu0 - A_puT)

				if dist_voisin>1.5*critere:
					ecartbis -= ecartbis/2
					precision_vectin -= precision_vectin/2.
				elif 0.5*critere>dist_voisin:
					ecartbis += ecartbis/2
					precision_vectin += precision_vectin/2.
				print "Point refstable : ", A[0,:]
				print "Point etustable : ", A[i,:]
				print "Point refinstable : ", B[0,:]
				print "Point etuinstable : ", Bbis
				print "Ecart : ", precision_vectin
				print "Distance entre les points : ", dist_voisin
				print "Critère : ", critere
				print "pour le point : ", i
				# time.sleep(5)

			# if dist_voisin<5*critere:
			results = np.vstack((results,A[i,:]))

	# Main loop :
	# -----------
	# for k in range(0,pt):
	# 	# Generating random variation around given direction :
	# 	# ----------------------------------------------------
	# 	# theta1 = rand.uniform(theta - espt, theta + espt)
	# 	# phi1   = rand.uniform(phi - espp, phi + espp)
	# 	theta1 = np.random.normal(theta, espt)
	# 	phi1   = np.random.normal(phi, espp)

	# 	# Computing direction :
	# 	# ---------------------
	# 	B = Fdirection(theta1,phi1,P_boundaries, Q_boundaries, U_boundaries, O, U_nominal, P_nominal)
	# 	# Computing last stable point in direction :
	# 	# ------------------------------------------
	# 	A[k,:] = Fexploration(O, B, Eps, seqPath, dtaPath, echPath, U_nominal, P_nominal)


	# Normalising points :
	# --------------------
	for i in range(0,results.shape[0]):
		results[i,:] = Conversion_si_pu(results[i,:],P_boundaries, Q_boundaries, U_boundaries)

	print "Matrice des points : ", results
	for i in range(1,results.shape[0]):
		print "Distance :", np.linalg.norm(results[i,:] - results[0,:])
	# G = np.mean(A, axis=0)

	# Computing covariance matrix :
	# -----------------------------
	M = np.cov(results.T)

	# Computing normal vector :
	# -------------------------
	n1 = np.linalg.svd(M)[0][:,-1]
	n2 = fitPLaneLTSQ(results)
	# n3 = fitPlaneEigen(A)
	n4 = fitPlaneSVD(results)
	# n5 = fitPlaneSolve(A)
	# n6 = np.cross(A[1,:] - A[0,:],A[2,:] - A[0,:])
	# n6 = n6 / np.linalg.norm(n6)

	print "Vct norm select n1 : ", n1
	print "Vct norm select n2 : ", n2
	# print "Vct norm select n3 : ", n3
	print "Vct norm select n4 : ", n4
	# print "Vct norm select n5 : ", n5
	# print "Vct norm select n6 : ", n6

	return results[0,:], n1

def FdomainBoundaries(parametersPath):
	""" 
		Compute study boundary depending on the size of the generator :
		Inputs :
			- parametersPath : Path to the parameters.par file
		Output :
			- Domain boundaries
			- Pnom (nominal Power) and Unom (nominal Tension)
	"""
	#Creating local variables :
	#--------------------------
	P_boundaries = np.zeros(2)
	Q_boundaries = np.zeros(2)
	U_boundaries = np.zeros(2)

	reg_name  = re.compile(r'(?<=name=")([^">]+)')
	reg_value = re.compile(r'(?<=value=")([^">]+)')
	name_Pnom = 'M2S_PNom'
	name_Unom = 'M2S_UNomGenTfo'

	#Reading parameters file :
	#-------------------------
	with open(parametersPath) as parametersFile:
		lines = parametersFile.readlines()

	for i, line in enumerate(lines):
		match = reg_name.search(str(line))
		if (match != None):
			if (match.group() == name_Pnom):
				Pnomapp = float(reg_value.findall(line)[0])
			if (match.group() == name_Unom): 
				Unom = float(reg_value.findall(line)[0])

	#Creating study domain : 
	#-------------------------
	P_boundaries[:] = [-1.5*Pnomapp,1.5*Pnomapp]
	Q_boundaries[:] = [-1.5*Pnomapp,1.5*Pnomapp]
	U_boundaries[:] = [0.5*Unom,1.5*Unom]

	# Retreiving reference name for mahcine :
	# ---------------------------------------
	# with open(MachinesDictionnary) as Dicofile:
	# 	for line in Dicofile:
	# 		if MachineName == line.split(";")[1]:
	# 			MachineRefName = line.split(";")[0]

	return P_boundaries, Q_boundaries, U_boundaries, Unom, Pnomapp

def Fdirection(theta, phi, P_boundaries, Q_boundaries, U_boundaries, O, U_nominal, P_nominal):
	"""
		Compute boundary point in specified direction 
		Inputs :
			- theta : angle (in rad) on the P,Q plane 
			- phi   : angle (in rad) between U axis and direction
			Note : Spherical coordinates
			- P_boundaries : Active power max and min
			- Q_boundaries : Reactive power max and min
			- U_boundaries : Tension max and min
			- O : Origin stable point	
		Output :
			- B : point on domain boundary

		Used in :
		- FplanTangent
	"""

	# Creating local vars :
	# ---------------------
	tf   = m.tan(phi)
	tt   = 0
	rho  = 1
	O_pu = Conversion_si_pu(O, P_boundaries, Q_boundaries, U_boundaries)
	#Creating output :
	#-----------------
	B = np.zeros(3)

	# Moving Origin & normalizing :
	# -----------------------------
	Pmax = 1 - O_pu[0]
	Pmin = 0 - O_pu[0]
	Qmax = 1 - O_pu[1]
	Qmin = 0 - O_pu[1]
	Umax = 1 - O_pu[2]
	Umin = 0 - O_pu[2]

	#Computing direction :
	#---------------------
	x = rho*m.sin(phi)*m.cos(theta)
	y = rho*m.sin(phi)*m.sin(theta)
	z = rho*m.cos(phi)

	# Computing coefficient to reach study border :
	# ---------------------------------------------
	if 0<=x:
		k_P = Pmax/x
	else:
		k_P = abs(Pmin/x)

	if 0<=y:
		k_Q = Qmax/y
	else:
		k_Q = abs(Qmin/y)

	if 0<=z:
		k_U = Umax/z
	else:
		k_U = abs(Umin/z)

	# Choose mininum coefficient and computing point on study border :
	# ----------------------------------------------------------------
	k = min(k_P, k_Q, k_U)
	B[:] = [k*x,k*y,k*z]

	# Moving back Origin :
	# --------------------
	B = B + O_pu

	# Converting to si from pu :
	# --------------------------
	B = Conversion_pu_si(B, P_boundaries, Q_boundaries, U_boundaries)

	return B

def Fraffinement(domain_data, faces, boundary_data, P_boundaries, Q_boundaries, U_boundaries):
	""" 
		Check if corners need to be refined. Update gross domain borders 
		Inputs :
			- points_pu : Domain corners
			- faces     : Domain faces information
			- A : boundary points 
		Output : 
			- corners : Corners to be refined
			- origins : Stable point close to corner
			- borders : Updated gross borders 
			- directions : Direction to reach corner from origin

		Used in :
			- main
	"""

	# Local variables :
	# -----------------
	index    = []
	epsilon  = 0.001
	epsilon2 = 0.00001
	p_plane  = np.zeros(shape=(3))
	error    = 0.

	# Outputs :
	# ---------
	directions = np.empty(shape=(0,2))
	origins    = np.empty(shape=(0,3))
	corners    = np.empty(shape=(0,3))

	# Updating boundaries :
	# ---------------------

	for i in range(0,len(domain_data)):
		# Finds which face includes point :
		# ---------------------------------
		index    = []
		index_ns = []
		# print "Point teste : ", points_pu[i,:]
		for j in range(0,len(faces)):
			
			# print "Face teste  : ", faces[j].points
			# print "TEST : ", np.all(np.isclose(faces[j].points,points_pu[i,:]),1).any()
			if np.all(np.isclose(faces[j].points,domain_data[i].point_pu),1).any():
				# Storing index :
				# ---------------
				index_ns.append(j)

		# print "nombre de point trouvé : ", len(index)

		if len(index_ns)>=3:
			print "Liste des indices : ", index_ns

			# Computing closest boundary points to corner :
			# ---------------------------------------------
			closest_points = []
			for k in range(0,len(index_ns)):
				print k, len(index_ns)
				closest_points.append((index_ns[k], np.linalg.norm(domain_data[i].point_pu - boundary_data[index_ns[k]].point_pu)))
			closest_points = sorted(closest_points, key=lambda closest_points: closest_points[1])
			index.append(closest_points[0][0])
			index.append(closest_points[1][0])
			index.append(closest_points[2][0])

			# Creating plane from boundary points :
			# -------------------------------------
			n_plane = np.cross((boundary_data[index[1]].point_pu-boundary_data[index[0]].point_pu),(boundary_data[index[2]].point_pu-boundary_data[index[0]].point_pu))
			n_plane = n_plane/np.linalg.norm(n_plane)

			# Computing projection point : 
			# ----------------------------
			d_plane = abs(np.dot(n_plane,domain_data[i].point_pu) - np.dot(n_plane,boundary_data[index[0]].point_pu))
			p_plane = PointDomain(domain_data[i].point_pu - d_plane*n_plane, 'pu', P_boundaries,Q_boundaries,U_boundaries)
			if abs(np.dot(p_plane.point_pu-boundary_data[index[0]].point_pu,n_plane))<epsilon2:
				p_plane = PointDomain(domain_data[i].point_pu - d_plane*n_plane, 'pu', P_boundaries,Q_boundaries,U_boundaries)
			else:
				p_plane = PointDomain(domain_data[i].point_pu + d_plane*n_plane, 'pu', P_boundaries,Q_boundaries,U_boundaries)


			# Checking if projection is in triangle :
			# ---------------------------------------
			aire  = np.linalg.norm(np.cross(boundary_data[index[2]].point_pu - boundary_data[index[1]].point_pu
										  , boundary_data[index[2]].point_pu - boundary_data[index[0]].point_pu))/2
			alpha = np.linalg.norm(np.cross(boundary_data[index[2]].point_pu - p_plane.point_pu
										  , boundary_data[index[1]].point_pu - p_plane.point_pu))/(2*aire)
			beta  = np.linalg.norm(np.cross(boundary_data[index[2]].point_pu - p_plane.point_pu
										  , boundary_data[index[0]].point_pu - p_plane.point_pu))/(2*aire)
			gamma = 1 - alpha - beta

			if ~(0<=alpha<=1 and 0<=beta<=1 and 0<=gamma<=1):
				# Computing gravitycentre :
				# -------------------------
				a = np.linalg.norm(boundary_data[index[2]].point_pu - boundary_data[index[1]].point_pu)
				b =	np.linalg.norm(boundary_data[index[2]].point_pu - boundary_data[index[0]].point_pu)
				c = np.linalg.norm(boundary_data[index[0]].point_pu - boundary_data[index[1]].point_pu)

				p_plane = PointDomain(boundary_data[index[0]].point_pu*a/(a+b+c) 
									+ boundary_data[index[1]].point_pu*b/(a+b+c)
									+ boundary_data[index[2]].point_pu*c/(a+b+c)
									, "pu",P_boundaries,Q_boundaries,U_boundaries)				

			d_plane = np.linalg.norm((domain_data[i].point_pu-p_plane.point_pu))
			print "Coin : ", domain_data[i].point_pu
			print "Proj : ", p_plane.point_pu
			print "La distance du coin au plan est : ", d_plane
			time.sleep(3)
			if error < d_plane:
				error = d_plane
			if d_plane>epsilon:

				# Computing direction :
				# ---------------------
				if (domain_data[i].point_pu[1]-p_plane.point_pu[1])>=0:
					theta = m.acos((domain_data[i].point_pu[0]-p_plane.point_pu[0]) 
									/ sqrt((domain_data[i].point_pu[0]-p_plane.point_pu[0])**2
										 + (domain_data[i].point_pu[1]-p_plane.point_pu[1])**2))
				else:
					theta = 2*pi - m.acos((domain_data[i].point_pu[0]-p_plane.point_pu[0]) 
									/ sqrt((domain_data[i].point_pu[0]-p_plane.point_pu[0])**2
										 + (domain_data[i].point_pu[1]-p_plane.point_pu[1])**2))
				# theta = m.atan((domain_data[i].point_pu[1]-p_plane.point_pu[1])/(domain_data[i].point_pu[0]-p_plane.point_pu[0]))
				# if theta<0:
				# 	theta = pi/2 - theta
				phi   = m.acos((domain_data[i].point_pu[2]-p_plane.point_pu[2])/np.linalg.norm(p_plane.point_pu - domain_data[i].point_pu))
				# Storing outputs : 
				# -----------------
				directions = np.vstack((directions, [theta,phi]))
				origins    = np.vstack((origins, p_plane.point_si))
				corners    = np.vstack((corners, domain_data[i].point_si))

	# Pmin = np.amin(corners[:,0])
	# Pmax = np.amax(corners[:,0])
	# Qmin = np.amin(corners[:,1])
	# Qmax = np.amax(corners[:,1])
	# Umin = np.amin(corners[:,2])
	# Umax = np.amax(corners[:,2])	
	# for i in range(0,len(boundary_data)):
	# 	Pmin = min(boundary_data[i].point_si[0],Pmin)
	# 	Pmax = max(boundary_data[i].point_si[0],Pmax)
	# 	Qmin = min(boundary_data[i].point_si[1],Qmin)
	# 	Qmax = max(boundary_data[i].point_si[1],Qmax)
	# 	Umin = min(boundary_data[i].point_si[2],Umin)
	# 	Umax = max(boundary_data[i].point_si[2],Umax)	

	return directions, origins, corners, error

def FcleanPoints(points, normals, epsilon):
	"""
		Clean the boundary points computed
		Inputs :
			- normals : matrix of normal vectors
			- points  : matrix of origin points
		Ouput  :
			- points_clean  : cleaned points
			- normals_clean : cleaned normals

		Used in :
			- main
	"""
	# Creating output :
	# -----------------
	points_clean  = np.zeros(shape=(0,3)) 
	normals_clean = np.zeros(shape=(0,3))

	# Local variable :
	# ----------------
	doublon = False

	for i in range(0,len(points)):
		for j in range(0,len(points)):
			diff_point = np.linalg.norm(points[i,:]-points[j,:])
			# Checking points position
			if (diff_point < epsilon):
				diff_normal = np.linalg.norm(np.cross(normals[i,:],normals[j,:]))
				# Checking points normal
				if (diff_normal < epsilon):
					doublon = True
		if ~doublon :
			points_clean  = np.vstack((points_clean,points[i,:]))
			normals_clean = np.vstack((normals_clean,normals[i,:]))

	return points_clean, normals_clean

# Ref domain :
# ------------
def FMonteCarlo(nbre_point, parametersPath, U_nominal, P_nominal ,P_boundaries,Q_boundaries,U_boundaries):
	"""
		Check for specified points if they are stable or not
		Inputs :
			- nbre_point : Number of points
		Ouput  :
			- points : the stable points
		Used in :
			- main
	"""
	# Creating output :
	# -----------------
	points_stable      = np.zeros(shape=(0,3))
	points_instable    = np.zeros(shape=(0,3))
	points_instable_in = np.zeros(shape=(0,3))
	points_unknown     = np.zeros(shape=(0,3))
	points_unknown_in  = np.zeros(shape=(0,3))
	outputsPath        = os.path.join(".","outputs","debug","globalInit","dumpInitValues-MACHINE.txt")

	# Main loop to check points stability :
	# -------------------------------------
	for i in range(0,nbre_point):
		# Computing random point in domain
		X_pu = [rand.uniform(0,1)
			, rand.uniform(0,1)
			, rand.uniform(0,1)]
		# print X_pu
		X = Conversion_pu_si(X_pu,P_boundaries,Q_boundaries,U_boundaries)

		# Checking Stability
		print "Point test : ", X
		X1 = FdynamoSimu(X ,parametersPath)
		alpha, var_instable = FreadOutputs(outputsPath)
		shutil.rmtree('outputs')
		print "Variable instable max : ", var_instable

		# If point is stable, it is added to output
		if (alpha == 0):
			points_stable = np.vstack((points_stable,X_pu))
		# elif (alpha == 1):
		# 	points_unknown = np.vstack((points_unknown,X_pu))
		elif (alpha >= 1):
			points_instable = np.vstack((points_instable,X_pu))

	hull = ConvexHull(points_stable)
	for i in range(0,points_instable.shape[0]):
		points_test = np.zeros(shape=(0,3))
		points_test = np.vstack((points_stable,points_instable[i,:]))
		hull_new = ConvexHull(points_test)
		if list(hull_new.vertices) == list(hull.vertices):
			points_instable_in = np.vstack((points_instable_in,points_instable[i,:]))

	# for i in range(0,points_unknown.shape[0]):
	# 	points_test = np.zeros(shape=(0,3))
	# 	points_test = np.vstack((points_stable,points_unknown[i,:]))
	# 	hull_new = ConvexHull(points_test)
	# 	if list(hull_new.vertices) == list(hull.vertices):
	# 		points_unknown_in = np.vstack((points_unknown_in,points_unknown[i,:]))

	return points_stable, points_instable_in # points_unknown_in

# Main Program :
# --------------
def main(version):

	# Parameters :
	# -----------
	p = 6 #number of points on boundary
	MachineRefName = "TEST_Ref"
	MachineName    = "TEST"

	# Path to data files :
	# --------------------
	parametersPath = os.path.join('parameters.par')

	#Creating data arrays :
	#----------------------
	boundary_data = []
	domain_data   = []
	A    = np.zeros(shape=(p,3)) # Store boundary points
	N    = np.zeros(shape=(p,3)) # Store boundary normal vector
	O    = np.zeros(shape=3) # Stable origin used
	O_pu = np.zeros(shape=3) # Stable origin in p.u 

	# Creating study domain from machine data :
	# -----------------------------------------
	P_boundaries, Q_boundaries, U_boundaries, U_nominal, P_nominal = FdomainBoundaries(parametersPath)
	
	print "Domaine d'etude :"
	print P_boundaries
	print Q_boundaries
	print U_boundaries, "\n"

	# Computing precision used in program from machine data :
	# -------------------------------------------------------
	print "Precision utilise :"
	pas = max(0.01,P_nominal*10**-4)
	precision_vect = Conversion_si_pu([P_boundaries[0] +  pas,Q_boundaries[0] + pas,U_boundaries[0] + pas*U_nominal/400.],P_boundaries, Q_boundaries, U_boundaries)
	print "Vecteur : ", precision_vect
	precision = np.linalg.norm(precision_vect)
	print "Norme : ", precision, "\n"

	# Choosing stable origin inside domaine : 
	# ---------------------------------------
	O[:] = [P_nominal*0.6,0,U_nominal]
	O_pu = Conversion_si_pu(O,P_boundaries,Q_boundaries,U_boundaries)

	# ======================================================= #
	#                         MAIN PROGRAM                    #
	# ======================================================= #

	# Version 0 : Computing stable simulation domain using tangent planes
	# ------------------------------------------------------------------- 
	if (version == 0):

		# Initialisation of domain :
		# --------------------------
		A[0,:], N[0,:] = FplanTangent(O, 0., pi/2, 9, P_boundaries, Q_boundaries, U_boundaries, precision, parametersPath, U_nominal, P_nominal, precision_vect)
		A[1,:], N[1,:] = FplanTangent(O, pi/2, pi/2, 9, P_boundaries, Q_boundaries, U_boundaries, precision, parametersPath, U_nominal, P_nominal, precision_vect)
		A[2,:], N[2,:] = FplanTangent(O, pi, pi/2, 9, P_boundaries, Q_boundaries, U_boundaries, precision, parametersPath, U_nominal, P_nominal, precision_vect)				
		A[3,:], N[3,:] = FplanTangent(O, 3*pi/2, pi/2, 9, P_boundaries, Q_boundaries, U_boundaries, precision, parametersPath, U_nominal, P_nominal, precision_vect)
		A[4,:], N[4,:] = FplanTangent(O, pi/4, 0.0, 9, P_boundaries, Q_boundaries, U_boundaries, precision, parametersPath, U_nominal, P_nominal, precision_vect)
		A[5,:], N[5,:] = FplanTangent(O, pi/4, pi, 9, P_boundaries, Q_boundaries, U_boundaries, precision, parametersPath, U_nominal, P_nominal, precision_vect)
		
		print "Points trouvés  : \n", A
		print "Vecteurs normaux : \n", N

		# Computing domain from support points and normals :
		faces, points = CreateMachineDomain.FdomainCreation(A,N,[0,1],[0,1],[0,1],O_pu)

		for i in range(0,p):
			boundary_data.append(PointBoundary(A[i,:],N[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))
		
		for i in range(0,len(points)):
			domain_data.append(PointDomain(points[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))

		# Creating export file
		# CreateMachineDomain.FdomainExport(faces, points, "test")
		# CreateMachineDomain.FconvexHullExport(A,"testconvHull")
		# CreateMachineDomain.FconvexHullExport(points,"testcornerHull")

		# Initialisation is done

		# Refinement step 1 :
		# -------------------
		directions, origins, corners, error = Fraffinement(domain_data, faces, boundary_data, P_boundaries, Q_boundaries, U_boundaries)

		# Computing new support points and normal for identified direction and origins :
		for i in range(0,len(directions)):
			print "Origine : ", origins[i,:]
			print "Directi : ", directions[i,:]
			print "Bord :", Fdirection(directions[i,0], directions[i,1], P_boundaries, Q_boundaries, U_boundaries, origins[i,:], U_nominal, P_nominal)
			new_pnt, new_vect = FplanTangent(origins[i,:], directions[i,0], directions[i,1], 9, P_boundaries, Q_boundaries, U_boundaries, precision, parametersPath, U_nominal, P_nominal, precision_vect)
			A = np.vstack((A,new_pnt))
			N = np.vstack((N,new_vect))

		# Computing domain from support points and normals :
		faces, points = CreateMachineDomain.FdomainCreation(A,N,[0,1],[0,1],[0,1],O_pu)

		for i in range(0,p):
			boundary_data.append(PointBoundary(A[i,:],N[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))
		
		for i in range(0,len(points)):
			domain_data.append(PointDomain(points[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))

		# Refinement step 2 :
		# -------------------
		directions, origins, corners, error = Fraffinement(domain_data, faces, boundary_data, P_boundaries, Q_boundaries, U_boundaries)

		# Computing new support points and normal for identified direction and origins :
		for i in range(0,len(directions)):
			print "Origine : ", origins[i,:]
			print "Directi : ", directions[i,:]
			print "Bord :", Fdirection(directions[i,0], directions[i,1], P_boundaries, Q_boundaries, U_boundaries, origins[i,:], U_nominal, P_nominal)
			new_pnt, new_vect = FplanTangent(origins[i,:], directions[i,0], directions[i,1], 9, P_boundaries, Q_boundaries, U_boundaries, precision, parametersPath, U_nominal, P_nominal, precision_vect)
			A = np.vstack((A,new_pnt))
			N = np.vstack((N,new_vect))

		# Computing domain from support points and normals :
		faces, points = CreateMachineDomain.FdomainCreation(A,N,[0,1],[0,1],[0,1],O_pu)

		for i in range(0,p):
			boundary_data.append(PointBoundary(A[i,:],N[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))
		
		for i in range(0,len(points)):
			domain_data.append(PointDomain(points[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))

	# Wrapping up simulation - compute error and create export :
	# ---------------------------------------------------------
	if (version == 0):
		# Converting output from internal units to (MW, Mvar, kV)
		corner_points_si   = np.zeros(shape=(points.shape[0],3))
		boundary_points_si = np.zeros(shape=(A.shape[0],3))
		for i in range(0,points.shape[0]):
			corner_points_si[i,:] = Conversion_pu_si(points[i,:], P_boundaries, Q_boundaries, U_boundaries)
		for i in range(0,A.shape[0]):
			boundary_points_si[i,:] = Conversion_pu_si(A[i,:], P_boundaries, Q_boundaries, U_boundaries)

		# Computing error between convex hull and constraints hull : 
		print "RAFFINEMENT POUR CALCUL DE LERREUR"
		directions, origins, corners, error = Fraffinement(domain_data, faces, boundary_data, P_boundaries,Q_boundaries,U_boundaries)
		print "FIN DU RAFFINEMENT POUR CALCUL DE LERREUR"

		print "Erreur Raffinement :", error

		# Creating export for convex hull in pu :
		CreateMachineDomain.FconvexHullExport(A,"testconvHull")
		CreateMachineDomain.FconvexHullExport(points,"testcornerHull")

		# Creating export for convex hull in si : 
		CreateMachineDomain.FconvexHullExport(boundary_points_si,"ConvHull_si")
		CreateMachineDomain.FexportCSV(boundary_points_si,"ampl_generators_domains",MachineRefName,MachineName,U_nominal,O)

		# Creating export for constraints hull : 
		CreateMachineDomain.FconvexHullExport(corner_points_si,"CornerHull_si")
		CreateMachineDomain.FexportCSV(corner_points_si,"ampl_generators_domains_coins",MachineRefName,MachineName,U_nominal,O)


if __name__ == '__main__':

	version = sys.argv[1]

	time_start = time.clock()
	main(int(version))
	time_elapsed = (time.clock() - time_start)
	print "Temps de calcul :", time_elapsed
	print "Fin du programme"
