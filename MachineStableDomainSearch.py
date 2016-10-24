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
from scipy.spatial import ConvexHull
import argparse

# Import Python visualisation libraries :
# ---------------------------------------
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib.path import Path



# Import custom dependencies :
# -----------------------------
import CreateMachineDomain

# Defining vector visualisation class :
# -------------------------------------
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

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
	B      = np.array([ [sum(X*Z), sum(Y*Z), sum(Z)] ])
	normal = np.linalg.solve(A,B.T)
	nn     = np.linalg.norm(normal)
	normal = normal / nn
	return normal.ravel()

def fitPlaneEigen(XYZ):
	average                  =sum(XYZ)/XYZ.shape[0]
	covariant                =np.cov(XYZ - average)
	eigenvalues,eigenvectors = np.linalg.eig(covariant)
	want_max                 = eigenvectors[:,eigenvalues.argmax()]
	(c,a,b)                  = want_max[3:6]
	normal                   = np.array([a,b,c])
	nn = np.linalg.norm(normal)
	return normal / nn  

# Program functions :
# -------------------		
def Fgroup(X, seqPath, dtaPath, echPath, P_boundaries, Q_boundaries, U_boundaries, U_nominal, P_nominal, precision):
	""" Call Eurostag simulator
		Inputs :
			- X       : machine parameters [P,Q,U]
			- seqPath : Machine .seq file 
			- dtaPath : Machine .dta file
			- echPath : Machine .ech file
			- Study boundaries
			- precision
		Output : 
			- alpha : Eurostag status
			- X1    : Eurostag stable point

		Used in :
			- Fexploration
			- Fraffinement
	"""

	# Local Variables : 
	# -----------------
	X1     = np.zeros(shape=3) # The simulation results are stored here
	E      = np.zeros(shape=3) # normalisation 
	t2     = 0.5 # Accepted biais
	alpha  = 0 # Initial indicator
	status = 5 # Eurostag status
	etat   = "ETAT D'EQUILIBRE A  0.10000D-02 VERIFIE POUR LES EQUATIONS MACHINE"
	etat2  = "ETAT D'EQUILIBRE A  0.10000D-02 NON VERIFIE DANS LES EQUATIONS SUIVANTES"

	# Setting up simulation :
	# -----------------------
	l.initLF(echPath)
	l.runLF()
	# print 'Load Flow succefull !'
	savPath = os.path.join("simTest.sav")
	l.initDynSimu(seqPath, dtaPath, savPath)

	# print 'Dynamic simulation initialised'
	MachineName = l.getMachines()
	Nodes = l.getNodes()
	# print "Machine :", MachineName, "\n"
	# print "Nodes :", Nodes, "\n"
	X1[2], angle = l.getNodeVoltage('N1      ', 'P')
	Unom, angle  = l.getNodeVoltage('N2      ', 'P')
	X1[0:2]      = l.getMachinePower(MachineName[0])
	print "P,Q,U INITIALEMENT :", X1

	# Checking Eurostag status :
	# --------------------------
	if etat in open("simTest.out").read():
		# Eurostag returns stable
		status = 0
		l.simulate()
	else:
		# Eurostag returns unstable :
		# ---------------------------
		equi_value = np.empty(0)
		with open("simTest.out") as outfile:
			lines = outfile.readlines()
			for i, line in enumerate(lines):
				if etat2 in line:
					k = i+5
					while (lines[k]!="\n"):
						value_line = lines[k].split(None)[-2]
						value_real = float(value_line.split("D")[0])*10**int(value_line.split("D")[1])
						equi_value = np.hstack((equi_value,value_real))
						k += 1
		# Checking tolerance :
		# --------------------
		# print equi_value
		if (np.amax(abs(equi_value)) < 0.01):
			status = 0
			l.simulate()
		elif (np.amax(abs(equi_value)) > 0.05):
			print "Cas problematique :"
			status = 1
			l.stopDynSimu()
		elif (0.01 < np.amax(abs(equi_value)) and np.amax(abs(equi_value)) < 0.5):
			status = 2
			l.simulate()
		else:
			print "Cas instable"
			status = 1
			l.simulate()
		# time.sleep(1)

	X1[2], angle = l.getNodeVoltage('N1      ', 'P')
	X1[0:2]      = l.getMachinePower(MachineName[0])
	print "P,Q,U FINALEMENT   :", X1	

	print ".... Eurostag status : ", status, " ....\n"

	print "#=====================================#"
	print "#    Fin de la simulation Eurostag    #"
	print "#=====================================#\n"
	# time.sleep(0.01)
	# Normalisation :
	# ---------------
	X_pu  =	Conversion_si_pu(X,P_boundaries,Q_boundaries,U_boundaries)
	X1_pu = Conversion_si_pu(X1,P_boundaries,Q_boundaries,U_boundaries)
	E = X_pu - X1_pu

	# Returning eurostag status :
	# ---------------------------
	print "Erreur :", np.linalg.norm(E)
	print "Preccs :", precision
	# time.sleep(1)
	if (np.linalg.norm(E)<precision*1000):
		if (status == 0):
			alpha = 0
		elif (status == 2):
			alpha = 1
		elif (status == 1):
			alpha = 2 
	else:
		alpha = 2

	return alpha, X1

def Fexploration(O,B,Esp,seqPath, dtaPath, echPath, U_nominal, P_nominal, P_boundaries, Q_boundaries, U_boundaries, MachineName):
	""" 
		Compute last stable point in a direction :
		Inputs : 
			- O   : Stable origin
			- B   : Boundary point in direction
			- Esp : Boundary precision
			- Study boundaries
			- Machine nominal power and tension
		Output :
			- A   : Last Stable point 

		Used in :
			- FplanTangent
	"""
	# Initialisation des variables internes :
	# ---------------------------------------
	C = np.zeros(3)          # point in unstable part
	A = np.zeros(3)          # point in stable part
	T = np.zeros(3)          # test point
	X = np.zeros(3)          # eurostag point result
	d = 0                    # distance between stable and unstable point
	# Initialisation :
	# ----------------
	A = O
	C = B
	T = C
	# distance has to be computed in pu !
	A_pu = Conversion_si_pu(A,P_boundaries,Q_boundaries,U_boundaries)
	C_pu = Conversion_si_pu(C,P_boundaries,Q_boundaries,U_boundaries)
	d = np.linalg.norm(A_pu-C_pu)

	# Last stable point at Eps biais :
	# --------------------------------
	while (d > Esp):
		print "Distance point stable/instable : ", d, "\n"

		it = 0
		# writing entry file for eurostag
		while (True and it<4):
			try:
				FeditechFile(T,"simTest",U_nominal, MachineName)
				break
			except IOError:
				it += 1
				time.sleep(0.5)
				pass

		print "Point test :", T
		try:
			[alpha,X] = Fgroup(T,seqPath, dtaPath, echPath, P_boundaries, Q_boundaries, U_boundaries, U_nominal, P_nominal, Esp)
		except:
			alpha = 3
			pass
		
		print "alpha =", alpha

		if (alpha == 0):
			#point is on stable side
			A = T             # New stable point
			T = (C+A)/2       # Updating test point
			print ".............. Cote Stable .............\n"
		else:
			#point is on unstable side
			C = T             # New unstable point
			T = (C+A)/2       # Updating test point
			print ".............. Cote instable .............\n"

		A_pu = Conversion_si_pu(A,P_boundaries,Q_boundaries,U_boundaries)
		C_pu = Conversion_si_pu(C,P_boundaries,Q_boundaries,U_boundaries)
		D_pu = C_pu - A_pu
		d = np.linalg.norm(D_pu)

		# Rest time for CPU :
		time.sleep(0.00000001)

	print "Dernier point stable   : ", A, '\n'
	print "Dernier point instable : ", C, '\n'
	return A

def Fvoisinage(B,P_boundaries,Q_boundaries, U_boundaries, sign1, sign2, precision_vect):
	""" 
		Return points in the vicinity of B
		Inputs :
			- B : point from which the vinicity is computed
			- Boundaries of domain
			- sign1/2 : + or -
		Output : 
			- voisins : Point in vicinity of B

		Used in :
			- FplanTangent
	"""	
	voisins    = np.zeros(shape=(3))  # Output
	voisins[:] = B                    # Initialisation
	aug_pas    = 100.                 # Bigger step

	[pas_P, pas_Q, pas_U] = Conversion_pu_si(precision_vect,P_boundaries,Q_boundaries,U_boundaries) - [P_boundaries[0], Q_boundaries[0], U_boundaries[0]]

	# B is on P boundary
	if np.any(np.isclose(B[0],P_boundaries)):
		voisins[1] += (1+sign1)*pas_Q*aug_pas
		voisins[2] += (1+sign2)*pas_U*aug_pas

	# B is on Q boundary
	elif np.any(np.isclose(B[1],Q_boundaries)):
		voisins[0] += (1+sign1)*pas_P*aug_pas
		voisins[2] += (1+sign2)*pas_U*aug_pas

	# B is on U boundary
	elif np.any(np.isclose(B[2],U_boundaries)):
		voisins[0] += (1+sign1)*pas_P*aug_pas
		voisins[1] += (1+sign2)*pas_Q*aug_pas

	# B is not on boundary
	else:
		return "Le point ", B, " n'est pas sur le bord du domaine"

	return voisins

def FplanTangent(O, theta, phi, pt, P_boundaries, Q_boundaries, U_boundaries, Eps, seqPath, dtaPath, echPath, U_nominal, P_nominal, precision_vect, MachineName):
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
	A = np.zeros(shape=(pt,3))        # Storage matrix of stable point
	B = np.zeros(shape=(pt,3))        # Storage matrix of unstable point
	results = np.zeros(shape=(0,3))   # Results matrix 
	M = np.zeros(shape=(3,3))         # Centered point matrix
	ecart = Eps*10                    # Step to compute other points in vicinity

	# method used to compute vicinity point :
	# 1 -> angle method (less stable)
	# 2 -> boundary step method (more stable)
	method = 2

	# Computing first point on stable domaine boundary :
	# --------------------------------------------------
	B[0,:] = Fdirection(theta,phi,P_boundaries, Q_boundaries, U_boundaries, O, U_nominal, P_nominal)
	A[0,:] = Fexploration(O, B[0,:], Eps, seqPath, dtaPath, echPath, U_nominal, P_nominal, P_boundaries, Q_boundaries, U_boundaries, MachineName)
	
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
				A[i,:]  = Fexploration(O, B[i,:], Eps, seqPath, dtaPath, echPath, U_nominal, P_nominal, P_boundaries, Q_boundaries, U_boundaries, MachineName)
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

			if dist_voisin<5*ecart:
				results = np.vstack((results,A[i,:]))

	if method == 2:
		precision_vectin = np.zeros(shape=(3))    # Store precision vector        
		critere  = Eps*100                       # Criteria for vicinity point
		ecartbis = critere                        # Step used to computed vicinity point
		Obis = O + 0.90*(O-A[0,:])                # New origin for faster exploration
		
		for i in range(1,9):
			precision_vectin[:] = precision_vect
			Bbis  = np.zeros(shape=(1,3))        # Point on boundary domain in the vicinity of B
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
			
			while (( (0.5*critere>dist_voisin or dist_voisin>1.5*critere) and it<=10) or it<1):
				it += 1
				Bbis   = Fvoisinage(B[0,:],P_boundaries,Q_boundaries,U_boundaries,sign1,sign2,precision_vectin)
				A[i,:] = Fexploration(O, Bbis, Eps, seqPath, dtaPath, echPath, U_nominal, P_nominal, P_boundaries, Q_boundaries, U_boundaries, MachineName)

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

			# if dist_voisin<5*critere:
			results = np.vstack((results,A[i,:]))

	# Normalising points :
	# --------------------
	for i in range(0,results.shape[0]):
		results[i,:] = Conversion_si_pu(results[i,:],P_boundaries, Q_boundaries, U_boundaries)

	print "Matrice des points : ", results
	for i in range(1,results.shape[0]):
		print "Distance :", np.linalg.norm(results[i,:] - results[0,:])

	# Computing covariance matrix :
	# -----------------------------
	# M = np.cov(results.T)
	M = np.vstack((results[:,0] - np.mean(results,axis=0)[0],
			results[:,1] - np.mean(results,axis=0)[1],
			results[:,2] - np.mean(results,axis=0)[2]))

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
	time.sleep(0.0000001)

	return results[0,:], n1

def FeditechFile(X,file,U_nominal, MachineName):
	""" 
		Create .ech Eurostag file (input Point of simulation)
		Inputs :
			- FileName : name of .ech File 
			- X = [P(MW), Q(Mvar), U(kV)] input point
		Output :
			- .ech file

		Used in :
			- Fexploration
	"""
	#Converting data to string :
	#---------------------------
	Pin  = str(X[0])[:8]
	Qin  = str(X[1])[:8]
	Uin  = str(X[2])[:8]
	Unom = str(U_nominal)[:8]
	
	#Setting up file path :
	#----------------------
	filename = file + ".ech"
	path     = os.path.join("")
	filename = os.path.join(path, filename)

	with open(filename, 'w') as echfile:
		#Writing .ech file :
		#------------------
		echfile.write("HEADER     "+ time.strftime("%d/%m/%y") +" 5.1\n")
		echfile.write(" \nB \n \n"
			"9 1 0 0 1 1 20    0.005                                 4   0         100. \n \n"
			"AA A \nAA B \n \n")
		echfile.write("1A N1                                                                             "
			+ Unom.rjust(10) +
			"            1.       0.       0.       0.\n"
			"1B N2                                                                             "
			+ Unom.rjust(10) +
			"            1.       0.       0.       0.\n")
		echfile.write(" \n5  N1                                        0.\n \n"
			"6 N1       N2      1       0.       0.\n \n"
			"G  " + MachineName.rjust(8) + " Y N1        -99999.       0.   99999.  -99999.       0.   99999. V"
			+ Uin.rjust(9) +
			" N1             1.       0.       0.\n \n")
		echfile.write("CH CHARGE   Y N2             0.       0."
			+ Pin.rjust(9) +
			"       0.       0."
			+ Qin.rjust(9) +
			"       0.       0.\n \n")

def FwritedtaFiles(MachinesdtaPath):
	""" 
		Create a specific dta file of each generators present in MachinesdtaPath
		*Note : not used (see IndusMachineDomaine.py) 
		Inputs  :
			- MachinesdtaPath : location of reference file
		Outputs : 
			- n single machine dta files

	"""

	path = os.path.join("dtaFiles")

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
				filename = MachineName + ".dta"
				filename = os.path.join(path, filename)
				with open(filename, 'w') as f:
					#Writing Header :
					#----------------
					f.write("HEADER     "+ time.strftime("%d/%m/%y") +" 5.1\n \n")
					#Writing Machine information :
					#-----------------------------
					while lines[i]!="\n":
						f.write(lines[i])
						i += 1
					#Writing Regulators information :
					#--------------------------------
					f.write(" \n")
					while (("M2       S") not in lines[i] and ("M2       U") not in lines[i]) and i<len(lines)-1:
						if lines[i].startswith("R " + MachineName):
							f.write(lines[i])
							f.write(lines[i+1])
							f.write(" \n \n \n")
						i += 1
					#Writing Network information :
					#-----------------------------
					f.write("I1 \nN2             0.     0.02      90.     100. \n \n \nLOADP   1\n         1              1.       1. \n \n \n \n \n"
						"CH \n1        W\n \n \n")

def FdomainBoundaries(dtaPath,MachinesDictionnary):
	""" 
		Extract generators caracteristics and create the study field from caracteristics
		Inputs :
			- dtaPath : location of generator parameters
			- MachinesDictionnary : location of dictionnary for generators
		Output :
			- Study boundaries (MW, Mvar, kV)
			- Unom : nominal tension (kV)
			- Pnomapp : Apparent nominal Power (MVar)
			- MachineName : Name of generator in Eurostag
			- MachineRefName : Name of generator in dictionnary file
		Used in : 
			- main
	"""
	#Creating local variables :
	#--------------------------
	P_boundaries = np.zeros(2)
	Q_boundaries = np.zeros(2)
	U_boundaries = np.zeros(2)
	MachineRefName = "None"

	#Opening specific machine dta file :
	#-----------------------------------
	with open(dtaPath) as Machinefile:
		#Reading file data :
		#-------------------
		lines    = Machinefile.readlines()
		#print(lines[3].split())
		# print(lines[3].split(" ")[4])
		MachineName = lines[3].split()[0]
		# Pnomturb = float(lines[6].split()[1])
		# Pnomalt  = float(lines[6].split()[2])
		
		try:
			Pnomapp = float(lines[3].split()[2])
			Unom    = float(lines[3].split()[3])
		except ValueError:
			Pnomapp = float(lines[3].split()[3])
			Unom    = float(lines[3].split()[4])
	#Creating machine domain : 
	#-------------------------
	P_boundaries[:] = [-1.5*Pnomapp,1.5*Pnomapp]
	Q_boundaries[:] = [-1.5*Pnomapp,1.5*Pnomapp]
	U_boundaries[:] = [0.5*Unom,1.5*Unom]

	# Finding MachineRefName :
	# ------------------------
	with open(MachinesDictionnary) as Dicofile:
		for line in Dicofile:
			# print MachineName
			if MachineName == line.split(";")[1]:
				MachineRefName = line.split(";")[0]

	return P_boundaries, Q_boundaries, U_boundaries, Unom, Pnomapp, MachineName, MachineRefName

def Fdirection(theta, phi, P_boundaries, Q_boundaries, U_boundaries, O, U_nominal, P_nominal):
	"""
		Compute boundary point in specified direction 
		Inputs :
			- theta : angle (in rad) on the P,Q plane 
			- phi   : angle (in rad) between U axis and direction
			Note : Spherical coordinates
			- P_boundaries, Q_boundaries, U_boundaries (MW, Mvar, kV) : Study field
			- O : Origin stable point of study
		Output :
			- B (MW, Mvar, kV) : point on domain boundary
		Used in :
		- FplanTangent
	"""

	# Creating local vars :
	# ---------------------
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
	x = m.sin(phi)*m.cos(theta)
	y = m.sin(phi)*m.sin(theta)
	z = m.cos(phi)

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

	k = min(k_P, k_Q, k_U)

	B[:] = [k*x,k*y,k*z]

	# Moving back Origin :
	# --------------------
	B = B + O_pu

	# Converting to si from pu :
	# --------------------------
	B = Conversion_pu_si(B, P_boundaries, Q_boundaries, U_boundaries)

	return B

def Fraffinement(domain_data,faces,boundary_data,P_boundaries,Q_boundaries,U_boundaries):
	""" 
		Check if corners need to be refined by projecting them on the convexHull of support points.
		The max distance of corner points to the convexHull of support point is also computed.
		Inputs :
			- domain_date : [n] list of PointDomain class containing corner points
			- faces : [n] list of faces class containing information of each faces 
			- boundary_data : [n] list of PointBoundary class containning support points
			- rest : not used
		Output : 
			- corners : Corners to be refined
			- origins : projection of corner on ConvexHull
			- directions : Direction to reach corner from origin
			- error_max : Max distance of corner points to ConvexHull

		Used in :
			- main (version 1 and 3)
			
		TO DO :
			- Correct upper bound of error between Convex hull and tangent planes domain (var : error_max)
	"""

	# Local variables :
	# -----------------
	epsilon   = 0.0000001             # criteria for belonging in plane
	p_plane   = np.zeros(shape=(3))   # point in plane
	error     = 0.                    
	error_max = 0.

	# Outputs :
	# ---------
	directions = np.empty(shape=(0,2))  # directions in theta and phi
	origins    = np.empty(shape=(0,3))  # origins for the direction
	corners    = np.empty(shape=(0,3))  # corresponding corner in direction

	for i in range(0,len(domain_data)):
		# Finds which face includes point :
		# ---------------------------------
		index    = []   # index of three closest support point
		index_ns = []   # index of all faces including point
		for j in range(0,len(faces)):
			
			if np.all(np.isclose(faces[j].points,domain_data[i].point_pu),1).any():
				index_ns.append(j)      # Storing index

		index_ns = np.arange(0,len(boundary_data))
		if len(index_ns)>=3:

			# Computing closest boundary points to corner :
			# ---------------------------------------------
			closest_points = []
			for k in range(0,len(index_ns)):
				closest_points.append([index_ns[k], np.linalg.norm(domain_data[i].point_pu - boundary_data[index_ns[k]].point_pu)])
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

			# finding real projection
			if abs(np.dot(p_plane.point_pu-boundary_data[index[0]].point_pu,n_plane))<epsilon:
				p_plane = PointDomain(domain_data[i].point_pu - d_plane*n_plane, 'pu', P_boundaries,Q_boundaries,U_boundaries)
			else:
				p_plane = PointDomain(domain_data[i].point_pu + d_plane*n_plane, 'pu', P_boundaries,Q_boundaries,U_boundaries)

			d_plane = np.linalg.norm((domain_data[i].point_pu-p_plane.point_pu))
			error   = d_plane

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
			print "Erreur first : ", error

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

			error = FcomputeError(domain_data[i].point_pu,boundary_data,error)
			if (error_max < error):
				error_max = error
			print "Erreur check : ", error_max

	return directions, origins, corners, error_max

def FquickDomainVisualisation(A, N, P_boundaries, Q_boundaries, U_boundaries,O,corners, origins):
	"""
		Create a cloud points to visualize next steps of computation :
		* Support points with normal vector
		* Corner points with their projection on the convexHull of the support points
		Inputs :
			- support_point : [n] list of PointBoundary Class
			- N : not used
			- P_boundaries, Q_boundaries, U_boundaries : boundaries of study (not used)
			- origin : stable origin used in the study
		Ouput  :
			- Visualisation of domain
		Used in :
			- main
	"""
	#Setting up visualiation :
	#-------------------------
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.set_aspect("equal")

	#Drawing domain boundaries :
	#---------------------------
	for s, e in combinations(np.array(list(product(P_boundaries,Q_boundaries,U_boundaries))), 2):
	    if np.sum(np.abs(s-e)) == (P_boundaries[1]-P_boundaries[0]):
	    	ax.plot3D(*zip(s,e), color="b")
	    elif np.sum(np.abs(s-e)) == (Q_boundaries[1]-Q_boundaries[0]):
	    	ax.plot3D(*zip(s,e), color="b")
	    elif np.sum(np.abs(s-e)) == (U_boundaries[1]-U_boundaries[0]):
	    	ax.plot3D(*zip(s,e), color="b")

	# Drawing vector associated to reference points :
	# -----------------------------------------------
	# for k in range(0,len(A)):
	# 	n = Arrow3D([A[k].point_si[0],A[k].point_si[0]+5*N[k,0]],[A[k].point_si[1],A[k].point_si[1]+5*N[k,1]],[A[k].point_si[2],A[k].point_si[2]+N[k,2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="r")	
	# 	ax.add_artist(n)

	# Drawing points :
	# ----------------
	colors = iter(cm.rainbow(np.linspace(0, 1, len(corners))))

	for k in range(0,len(A)):
		ax.scatter(A[k].point_si[0],A[k].point_si[1],A[k].point_si[2], color='g',s=20)
	for k in range(0,len(corners)):
		c = next(colors)
		ax.scatter(corners[k,0],corners[k,1],corners[k,2],color=c, zdir='A[:,2]',s=50)
		ax.scatter(origins[k,0],origins[k,1],origins[k,2],color=c, zdir='A[:,2]',s=50)
		a = Arrow3D([origins[k,0],corners[k,0]],[origins[k,1],corners[k,1]],[origins[k,2],corners[k,2]], mutation_scale=10, lw=1, arrowstyle="-|>", color=c)
		ax.add_artist(a)
		plt.draw()

	ax.scatter(O[0],O[1],O[2],color='g',s=2)
	plt.draw()
	#Display figure :
	#----------------
	ax.set_xlabel('P', fontsize='xx-large')
	ax.set_ylabel('Q', fontsize='xx-large')
	ax.set_zlabel('U', fontsize='xx-large')
	# fig.suptitle(MachineName, fontsize='xx-large')
	fig.savefig('QuickViz')
	plt.draw()

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

def FcomputeError(corner_point,boundary_data,error):
	"""
		Check if previous calculated error is the minimum distance to boundary domain
		Inputs :
			- corner_point : studied point
			- boundary_data : points of the boundary
			- error : distance between corner_point and its custom projection on boundary domain
		Outputs :
			dist_point : minimum distance of corner_point to the boundary domain
	"""

	dist_point = error   # Initialising dist_point

	# Checking the distance of corner_point from each boundary points
	for i in range(0,len(boundary_data)):
		if np.linalg.norm(corner_point - boundary_data[i].point_pu)<dist_point:
			# Updating dist_point if a new min is found
			dist_point = np.linalg.norm(corner_point - boundary_data[i].point_pu)

	return dist_point 

def FMonteCarlo(nbre_point,seqPath, dtaPath, echPath, U_nominal, P_nominal ,P_boundaries,Q_boundaries,U_boundaries, precision, MachineName):
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
	points_stable   = np.zeros(shape=(0,3))
	points_instable = np.zeros(shape=(0,3))
	points_instable_in = np.zeros(shape=(0,3))
	points_unknown = np.zeros(shape=(0,3))
	points_unknown_in = np.zeros(shape=(0,3))

	# Main loop to check points stability :
	# -------------------------------------
	for i in range(0,nbre_point):
		# Computing random point in domain
		X_pu = [rand.uniform(0.,1.)
			, rand.uniform(0.,1.)
			, rand.uniform(0.,1.)]
		# print X_pu
		X = Conversion_pu_si(X_pu,P_boundaries,Q_boundaries,U_boundaries)
		# Checking Stability
		FeditechFile(X,"simTest",U_nominal, MachineName)

		print "Point test : ", X
		try:
			alpha, X1 = Fgroup(X ,seqPath, dtaPath, echPath, P_boundaries, Q_boundaries, U_boundaries, U_nominal, P_nominal, precision)
		except:
			alpha = 3
			pass
		# FeditechFile(X,str(X) + "_" + str(alpha),U_nominal)
		# If point is stable, it is added to output
		if (alpha == 0):
			points_stable = np.vstack((points_stable,X_pu))
		elif (alpha == 1):
			points_unknown = np.vstack((points_unknown,X_pu))
		elif (alpha > 1):
			points_instable = np.vstack((points_instable,X_pu))

	hull = ConvexHull(points_stable)
	for i in range(0,points_instable.shape[0]):
		points_test = np.zeros(shape=(0,3))
		points_test = np.vstack((points_stable,points_instable[i,:]))
		hull_new = ConvexHull(points_test)
		if list(hull_new.vertices) == list(hull.vertices):
			points_instable_in = np.vstack((points_instable_in,points_instable[i,:]))
	for i in range(0,points_unknown.shape[0]):
		points_test = np.zeros(shape=(0,3))
		points_test = np.vstack((points_stable,points_unknown[i,:]))
		hull_new = ConvexHull(points_test)
		if list(hull_new.vertices) == list(hull.vertices):
			points_unknown_in = np.vstack((points_unknown_in,points_unknown[i,:]))

	return points_stable, points_instable_in, points_unknown_in

# Main Program :
# --------------
def main(version, visualisation):
	"""
	Main program, compute 3D stable domain of generator :
	Inputs : 
		- version : which version of the program is used 
		0 : simple   (Initialisation + 2*refinement) 
		1 : sampling (random points)
		2 : write specified ech file for manual simulation
		3 : criteria based on error and max iteration
		4 : simple+  (Initiation + 2*refinement + corners)
		- visualisation : on or off
		0 : off
		1 : on
	Ouputs :
		- for each stage a .off output is created
		- a .csv file is created containing all planes
	Used in : 
		- IndusMachineDomaine.py
	"""

	# Setting up visualization :
	# --------------------------
	if (visualisation == 1 or visualisation == 0):
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

	# Parameters :
	# -----------
	p = 6               # number of points on boundary
	iteration = 0       # number of iteration
	error = 0           # initialisation of error

	# Path to data files :
	# --------------------
	dtaFile = os.path.join("simTest.dta")
	echFile = os.path.join("simTest.ech")
	seqFile = os.path.join("simTest.seq")
	
	#Creating data arrays :
	#----------------------
	boundary_data = []
	domain_data   = []
	A    = np.zeros(shape=(p,3)) 	# Store boundary points
	N    = np.zeros(shape=(p,3)) 	# Store boundary normal vector
	O    = np.zeros(shape=3) 		# Stable origin used
	O_pu = np.zeros(shape=3) 		# Stable origin in p.u 

	# Creating study domain from machine data :
	# -----------------------------------------
	P_boundaries, Q_boundaries, U_boundaries, U_nominal, P_nominal, MachineName, MachineRefName = FdomainBoundaries(dtaFile, DictPath)
	
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
	if (version == 0 or version == 4):
		if (version == 0):	
			print "VERSION : Simple"
		if (version == 4):	
			print "VERSION : Simple+"

		print "#==================================#"
		print "#                                  #"
		print "#          INITIALISATION          #"
		print "#                                  #"
		print "#==================================#"

		A[0,:], N[0,:] = FplanTangent(O, 0., pi/2, 9, P_boundaries, Q_boundaries, U_boundaries, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, precision_vect, MachineName)
		A[1,:], N[1,:] = FplanTangent(O, pi/2, pi/2, 9, P_boundaries, Q_boundaries, U_boundaries, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, precision_vect, MachineName)
		A[2,:], N[2,:] = FplanTangent(O, pi, pi/2, 9, P_boundaries, Q_boundaries, U_boundaries, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, precision_vect, MachineName)				
		A[3,:], N[3,:] = FplanTangent(O, 3*pi/2, pi/2, 9, P_boundaries, Q_boundaries, U_boundaries, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, precision_vect, MachineName)
		A[4,:], N[4,:] = FplanTangent(O, pi/4, 0.0, 9, P_boundaries, Q_boundaries, U_boundaries, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, precision_vect, MachineName)
		A[5,:], N[5,:] = FplanTangent(O, pi/4, pi, 9, P_boundaries, Q_boundaries, U_boundaries, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, precision_vect, MachineName)
		

		print "Points trouvés  : \n", A
		print "Vecteurs normaux : \n", N

		print " \n INITIALISATION TERMINE \n"

		for i in range(0,p):
			boundary_data.append(PointBoundary(A[i,:],N[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))

		# Creating domain
		faces, points = CreateMachineDomain.FdomainCreation(A,N,[0,1],[0,1],[0,1],O_pu,precision)
		
		for i in range(0,len(points)):
			domain_data.append(PointDomain(points[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))

		# Export domain
		CreateMachineDomain.FconvexHullExport(A,"testconvHull_0")
		CreateMachineDomain.FconvexHullExport(points,"testcornerHull_0")
		
		print "#==================================#"
		print "#                                  #"
		print "#         1er RAFFINEMENT          #"
		print "#                                  #"
		print "#==================================#"

		directions, origins, corners, error = Fraffinement(domain_data, faces, boundary_data, P_boundaries,Q_boundaries,U_boundaries)

		# Rest time for CPU
		time.sleep(0.00001)

		if (visualisation == 1):
			FquickDomainVisualisation(boundary_data,N,P_boundaries,Q_boundaries,U_boundaries,O,corners,origins)

		for i in range(0,len(directions)):
			new_pnt, new_vect = FplanTangent(origins[i,:], directions[i,0], directions[i,1], 9, P_boundaries, Q_boundaries, U_boundaries, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, precision_vect, MachineName)
			A = np.vstack((A,new_pnt))
			N = np.vstack((N,new_vect))
			# Rest time for CPU
			time.sleep(0.00001)

		print " \n 1er RAFFINEMENT TERMINE \n"

		# Creating domain
		faces, points = CreateMachineDomain.FdomainCreation(A,N,[0,1],[0,1],[0,1],O_pu,precision)

		# Export domain
		CreateMachineDomain.FconvexHullExport(A,"testconvHull_1")
		CreateMachineDomain.FconvexHullExport(points,"testcornerHull_1")

		boundary_data = []
		for i in range(0,len(A)):
			boundary_data.append(PointBoundary(A[i,:],N[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))

		domain_data = []
		for i in range(0,len(points)):
			domain_data.append(PointDomain(points[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))

		print "#==================================#"
		print "#                                  #"
		print "#         2nd RAFFINEMENT          #"
		print "#                                  #"
		print "#==================================#"

		directions, origins, corners, error = Fraffinement(domain_data, faces, boundary_data, P_boundaries,Q_boundaries,U_boundaries)
		
		# Rest time for CPU
		time.sleep(0.00001)

		if (visualisation == 1):
			FquickDomainVisualisation(boundary_data,N,P_boundaries,Q_boundaries,U_boundaries,O,corners,origins)

		for i in range(0,len(directions)):

			new_pnt, new_vect = FplanTangent(origins[i,:], directions[i,0], directions[i,1], 9, P_boundaries, Q_boundaries, U_boundaries, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, precision_vect, MachineName)
			A = np.vstack((A,new_pnt))
			N = np.vstack((N,new_vect))
			# Rest time for CPU
			time.sleep(0.00001)

		print " \n 2nd RAFFINEMENT TERMINE \n"

		# Creating domain
		faces, points = CreateMachineDomain.FdomainCreation(A,N,[0,1],[0,1],[0,1],O_pu,precision)

		# Export domain
		CreateMachineDomain.FconvexHullExport(A,"testconvHull_2")
		CreateMachineDomain.FconvexHullExport(points,"testcornerHull_2")

		boundary_data = []
		for i in range(0,len(A)):
			boundary_data.append(PointBoundary(A[i,:],N[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))

		domain_data = []
		for i in range(0,len(points)):
			domain_data.append(PointDomain(points[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))

	# Version 1 : Computing stable simulation domain using sampling
	# the number of points can be specified in FMonteCarlo
	# -------------------------------------------------------------
	if (version == 1):
		print "VERSION : Aleatoire"

		points_MC_stable, points_MC_instable, points_MC_unknown = FMonteCarlo(15000,seqFile, dtaFile, echFile, U_nominal, P_nominal,P_boundaries,Q_boundaries,U_boundaries,precision, MachineName)
		CreateMachineDomain.FconvexHullExport(points_MC_stable,"testconvHull_MC")
		
		print "Fin des simulations"

		points_MC_stable_OK = points_MC_stable
		for i in range(0,points_MC_instable.shape[0]):
			j = 0
			pointins_inDomain = True
			while pointins_inDomain and j<points_MC_stable_OK.shape[0]:
				points_test = np.delete(points_MC_stable_OK, j, 0)
				hull = ConvexHull(points_test)
				points_test = np.vstack((points_test,points_MC_instable[i,:]))
				hull_new = ConvexHull(points_test)
				if list(hull_new.vertices) != list(hull.vertices):
					pointins_inDomain = False
					points_MC_stable_OK = np.delete(points_MC_stable_OK, j, 0) 
				j+=1

		points_MC_stable_si    = np.zeros(shape=(points_MC_stable.shape[0],3))
		points_MC_instable_si  = np.zeros(shape=(points_MC_instable.shape[0],3))
		points_MC_unknown_si   = np.zeros(shape=(points_MC_unknown.shape[0],3))
		points_MC_stable_OK_si = np.zeros(shape=(points_MC_stable_OK.shape[0],3))

		for i in range(0,points_MC_stable.shape[0]):
			points_MC_stable_si[i,:] = Conversion_pu_si(points_MC_stable[i,:], P_boundaries, Q_boundaries, U_boundaries)
		for i in range(0,points_MC_instable.shape[0]):
			points_MC_instable_si[i,:] = Conversion_pu_si(points_MC_instable[i,:], P_boundaries, Q_boundaries, U_boundaries)
			point_filename = MachineName + "_" +str(i)
			# FeditechFile(points_MC_instable_si[i,:],point_filename,U_nominal, MachineName)   # Activate this line to export non-convex unstable points
		for i in range(0,points_MC_unknown.shape[0]):
			points_MC_unknown_si[i,:] = Conversion_pu_si(points_MC_unknown[i,:], P_boundaries, Q_boundaries, U_boundaries)
		for i in range(0,points_MC_stable_OK.shape[0]):
			points_MC_stable_OK_si[i,:] = Conversion_pu_si(points_MC_stable_OK[i,:], P_boundaries, Q_boundaries, U_boundaries)			

		hull = ConvexHull(points_MC_stable_OK_si)

		# Export domain
		CreateMachineDomain.FexportCSV(points_MC_stable_OK_si,"ampl_generators_domains_MC",MachineRefName,MachineName,U_nominal,O)

		print "Points stables   : ", points_MC_stable.shape[0]
		print "Points instables : ", points_MC_instable.shape[0]
		print "Points inconnus  : ", points_MC_unknown.shape[0]
		print "Points stables ok: ", points_MC_stable_OK.shape[0]

		if (visualisation == 1 or visualisation == 0):
			plt.plot(points_MC_stable_si[:,0],points_MC_stable_si[:,1],points_MC_stable_si[:,2], 'o', color='orange')
			plt.plot(points_MC_instable_si[:,0],points_MC_instable_si[:,1],points_MC_instable_si[:,2], 'o', color="red")
			plt.plot(points_MC_unknown_si[:,0],points_MC_unknown_si[:,1],points_MC_unknown_si[:,2], 'o', color="yellow")
			plt.plot(points_MC_stable_OK_si[:,0],points_MC_stable_OK_si[:,1],points_MC_stable_OK_si[:,2], 'o', color='green')
			
			for simplex in hull.simplices:
				plt.plot(points_MC_stable_OK_si[simplex,0], points_MC_stable_OK_si[simplex,1], points_MC_stable_OK_si[simplex,2], 'k-')

			ax.set_xlabel('P', fontsize='xx-large')
			ax.set_ylabel('Q', fontsize='xx-large')
			ax.set_zlabel('U', fontsize='xx-large')
			fig.suptitle(MachineName, fontsize='xx-large')
			fig.savefig(MachineName + '_MC.png')

			plt.draw()

	# Version 2 : Used to export problematic .ech files :
	# ---------------------------------------------------
	if (version == 2):
		T1 = np.zeros(3)
		T2 = np.zeros(3)
		T3 = np.zeros(3)
		T4 = np.zeros(3)

		T1[:] = [726.59997559, 944.57995605, 24.]
		T2[:] = [726.59997559, -1198.89001465, 24.]
		T3[:] = [727.81097412, 1.21099997, 32.1600008]
		T4[:] = [726.59997559, 0., 13.20000029]

		FeditechFile(T1,"simTest_Qmax",U_nominal,MachineName)
		FeditechFile(T2,"simTest_Qmin",U_nominal,MachineName)
		FeditechFile(T3,"simTest_Umax",U_nominal,MachineName)
		FeditechFile(T4,"simTest_Umin",U_nominal,MachineName)

	# Version 3 : Computing stable simulation domain using tangent planes and adding 
	# a specific confition to stop refinement (5 step max)
	# ------------------------------------------------------------------------------
	if (version == 3):
		print "VERSION : Critere"
		A[0,:], N[0,:] = FplanTangent(O, 0., pi/2, 9, P_boundaries, Q_boundaries, U_boundaries, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, precision_vect, MachineName)
		A[1,:], N[1,:] = FplanTangent(O, pi/2, pi/2, 9, P_boundaries, Q_boundaries, U_boundaries, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, precision_vect, MachineName)
		A[2,:], N[2,:] = FplanTangent(O, pi, pi/2, 9, P_boundaries, Q_boundaries, U_boundaries, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, precision_vect, MachineName)				
		A[3,:], N[3,:] = FplanTangent(O, 3*pi/2, pi/2, 9, P_boundaries, Q_boundaries, U_boundaries, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, precision_vect, MachineName)
		A[4,:], N[4,:] = FplanTangent(O, pi/4, 0.0, 9, P_boundaries, Q_boundaries, U_boundaries, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, precision_vect, MachineName)
		A[5,:], N[5,:] = FplanTangent(O, pi/4, pi, 9, P_boundaries, Q_boundaries, U_boundaries, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, precision_vect, MachineName)

		for i in range(0,p):
			boundary_data.append(PointBoundary(A[i,:],N[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))

		faces, points = CreateMachineDomain.FdomainCreation(A,N,[0,1],[0,1],[0,1],O_pu,precision)
		
		for i in range(0,len(points)):
			domain_data.append(PointDomain(points[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))

		# Export Domain
		CreateMachineDomain.FconvexHullExport(A,"testconvHull")
		CreateMachineDomain.FconvexHullExport(points,"testcornerHull")

		while (error>0.05 or iteration<1) and (iteration<3):
			iteration += 1
			directions, origins, corners, error = Fraffinement(domain_data, faces, boundary_data, P_boundaries,Q_boundaries,U_boundaries)
			for i in range(0,len(directions)):
				new_pnt, new_vect = FplanTangent(origins[i,:], directions[i,0], directions[i,1], 9, P_boundaries, Q_boundaries, U_boundaries, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, precision_vect, MachineName)
				A = np.vstack((A,new_pnt))
				N = np.vstack((N,new_vect))

			# Creating domain
			faces, points = CreateMachineDomain.FdomainCreation(A,N,[0,1],[0,1],[0,1],O_pu,precision)

			# Export Domain
			CreateMachineDomain.FconvexHullExport(A,"testconvHull")
			CreateMachineDomain.FconvexHullExport(points,"testcornerHull")

			boundary_data = []
			for i in range(0,len(A)):
				boundary_data.append(PointBoundary(A[i,:],N[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))

			domain_data = []
			for i in range(0,len(points)):
				domain_data.append(PointDomain(points[i,:],"pu",P_boundaries,Q_boundaries,U_boundaries))

	# Preparing and writing finale results :
	# --------------------------------------
	if (version == 0 or version == 3 or version == 4):

		# Compute error
		directions, origins, corners, error = Fraffinement(domain_data, faces, boundary_data, P_boundaries,Q_boundaries,U_boundaries)

		# Computing new stables point in version simple +
		if (version == 4):
			for i in range(0,len(directions)):
				border_point = Fdirection(directions[i,0], directions[i,1],P_boundaries, Q_boundaries, U_boundaries, origins[i,:], U_nominal, P_nominal)
				new_pnt = Fexploration(origins[i,:], border_point, precision, seqFile, dtaFile, echFile, U_nominal, P_nominal, P_boundaries, Q_boundaries, U_boundaries, MachineName)
				A = np.vstack((A,Conversion_si_pu(new_pnt,P_boundaries, Q_boundaries, U_boundaries)))
				# Rest time for CPU
				time.sleep(0.00001)

		# Converting point from pu to si for export 
		corner_points_si   = np.zeros(shape=(points.shape[0],3))
		boundary_points_si = np.zeros(shape=(A.shape[0],3))
		for i in range(0,points.shape[0]):
			corner_points_si[i,:] = Conversion_pu_si(points[i,:], P_boundaries, Q_boundaries, U_boundaries)
		for i in range(0,A.shape[0]):
			boundary_points_si[i,:] = Conversion_pu_si(A[i,:], P_boundaries, Q_boundaries, U_boundaries)

		hull = ConvexHull(boundary_points_si)      # Convex hull for visualisation 

		if (visualisation == 1 or visualisation == 0):
			plt.plot(boundary_points_si[:,0],boundary_points_si[:,1],boundary_points_si[:,2], 'o')

			for simplex in hull.simplices:
				plt.plot(boundary_points_si[simplex,0], boundary_points_si[simplex,1], boundary_points_si[simplex,2], 'k-')

			ax.set_xlabel('P', fontsize='xx-large')
			ax.set_ylabel('Q', fontsize='xx-large')
			ax.set_zlabel('U', fontsize='xx-large')
			fig.suptitle(MachineName, fontsize='xx-large')
			fig.savefig(MachineName + '_MC.png')
			plt.draw()


		if (visualisation == 1):
			FquickDomainVisualisation(boundary_data,N,P_boundaries,Q_boundaries,U_boundaries,O,corners,origins)

		# Writing outputs and export domain and constraints :
		print "Erreur Raffinement :", error
		CreateMachineDomain.FconvexHullExport(A,"ConvHull_final")
		CreateMachineDomain.FexportCSV(boundary_points_si,"ampl_generators_domains",MachineRefName,MachineName,U_nominal,O)
		CreateMachineDomain.FconvexHullExport(points,"CornerHull_final")
		CreateMachineDomain.FexportCSV(corner_points_si,"ampl_generators_domains_coins",MachineRefName,MachineName,U_nominal,O)

if __name__ == '__main__':	

	# Import Eurostag DLL :
	# ---------------------
	eurostag_path = os.environ.get('EUROSTAG')
	if (eurostag_path is None):
		eurostag_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Eurostag"))

	try:
		l = imp.load_source ('eurostag', os.path.join(eurostag_path, 'python_interface', 'eurostag.py'))
	except Exception as exc:
		print("Failed to load Eurostag API :" + str(exc))
		sys.exit(1)

	print "Eurostag has been succesfully loaded\n"

	parser = argparse.ArgumentParser(description='Compute stable simulation domain for generator in folder. To run the program, you must be first be inside the simulation directory and run "python ..\MachineStableDomainSearch.py + arguments"')
	parser.add_argument('-v','--version', help='which version of the program to be run : - 0 Simple, - 1 Random, - 3 Criteria,  4 - Simple+',
			 required=True)
	parser.add_argument('-i','--image', help='switch visualization on or off (0 off, 1 on)',
			 required=True)
	args = parser.parse_args()


	DictPath   = os.path.join("..","dict.csv")

	time_start = time.clock()
	main(int(args.version), int(args.image))
	time_elapsed = (time.clock() - time_start)
	print "Temps de calcul :", time_elapsed
	print "Fin du programme"
	if (int(args.image) == 1):
		plt.show()

	l.finalize() 	
