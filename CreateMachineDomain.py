
# -*- coding: utf-8 -*-

# Copyright (c) 2016, Pierre Saikaly  (saikalypierre@gmail.com)
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#===========================#
# created on 20 april 2016
#===========================#

#Import Python dependencies :
#----------------------------
import os
import numpy as np
import random as rand
from numpy import linalg
import math as m
from math import pi
from math import sqrt
from math import fabs
import time

from scipy.spatial import ConvexHull
from scipy.optimize import curve_fit

# Import Python visualisation dependencies :
# -----------------------------------------
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class PointDomain():
	"""Define point of domain in pu and si"""

	def __init__(self,p_input,unit):
		"""
			Define point and normal in both pu and si :
				Inputs :
					p_input : P,Q,U point
					unit    : si or pu
				Output :
					p_boundary_si : P,Q,U in si
					p_boundary_pu : P,Q,U in pu
		"""
		if unit=="pu":
			self.p_boundary_pu = p_input
		elif unit=="si":
			self.p_boundary_si = p_input
		else:
			print "Specifier l'unite : pu ou si"

class Face:
	"""Define faces by adding points and ordering them for OFF format."""

	def __init__(self, normal):
		"""By default, the faces is just defined by it's normal vector"""
		self.points = np.empty(shape=(0,3))
		self.normal = normal

	def addPoint(self, point):
		self.points = np.vstack((self.points,point))

	def orderPoints(self):
		"""
			By default when the points are added to the surface they are not 
			ordered. This method orders the point for OFF format output
		"""
		# Local variable :
		# ----------------
		point_ordered = np.empty(shape=(0,3))
		point_ordered = np.vstack((point_ordered,self.points[0,:])) # starting point is first point
		self.points   = np.delete(self.points, 0, 0) 

		# Main loop for ordering points :
		# ------------------------------- 
		while len(self.points) != 0:
			last_point  = point_ordered[-1,:]
			first_point = self.points[0,:]
			for i in range(0,len(self.points)):
				v1 = first_point - last_point
				v2 = self.points[i,:] - last_point
				v3 = np.cross(v1,v2)
				if np.dot(self.normal,v3)<0:
					first_point = self.points[i,:]
			# Adding ordered point
			point_ordered = np.vstack((point_ordered,first_point))
			# Removing ordered point
			self.points = np.delete(self.points,np.where(np.all(self.points==first_point,axis=1)),0)



		# Storing back ordered points :
		# -----------------------------
		self.points = point_ordered

def FinDomain(point,A,N,O):
	"""
		Check if point is in domain
		Inputs :
			- point : Test point
			- N : matrix of normal vectors
			- A : matrix of origin points
			- O : stable origin
		Ouput  :
			- True or False
		Used in :
			- FdomainCreation
	"""	

	# Local Variable :
	# ----------------
	epsilon = 0.0000000001   # Tolerance
	# Testing point : 
	# ---------------
	for i in range(0,len(A)):
		if (np.dot( (point-A[i,:]),N[i,:]) * np.dot((O-A[i,:]),N[i,:]) / abs(np.dot((O-A[i,:]),N[i,:]) ) < -epsilon):
			return False
			break

	return True

def FchecknormalPlane(points_support,normals,origin,precision):
	"""
		Check if normals include all support points. If it does not, the normal is
		removed from the list defining the domain.
		Inputs :
			- normals : matrix of normal vectors
			- points_support : matrix of support points
			- origin : Stable origin
		Ouput  :
			- points_cleaned
			- normals_cleaned
		Used in :
			- FdomainCreation
	"""		
	epsilon = precision*100   # Tolerance 
	points_cleaned  = np.zeros(shape=(0,3)) 
	normals_cleaned = np.zeros(shape=(0,3))

	for j in range(0,len(normals)):
		normalcorrect = True
		for i in range(0,len(points_support)):
			count = 0
			# Checking if point is in domain
			if (np.dot( (points_support[i,:]-points_support[j,:]),normals[j,:]) * np.dot((origin-points_support[j,:]),normals[j,:]) 
				/ abs(np.dot((origin-points_support[j,:]),normals[j,:]) ) < -epsilon):
					normalcorrect = False # Normal is flagged as incorrect
					count += 1 # Updating number of points not in domain
			print "Nombre de vecteurs faux : ", count
		if normalcorrect :
			# Adding point and normal
			points_cleaned = np.vstack((points_cleaned,points_support[j,:]))
			normals_cleaned = np.vstack((normals_cleaned, normals[j,:]))



	return points_cleaned, normals_cleaned

def FdomainCreation(points_support,normals,P_boundaries, Q_boundaries, U_boundaries,O, precision):
	"""
		Create the faces and vertex of domain
		Inputs :
			- normal : matrix of normal vectors
			- points_support : matrix of support points
			- boundaries : study boundaries
			- O : stable origin
			- precision : 
		Ouput  :
			- faces_   : composed of normal vector and ordered vertex
			- points_ : contain all the vertex of the domain
		Used in :
			- 
	"""	

	# Cleaning points :
	# -----------------
	A, N = FchecknormalPlane(points_support,normals,O, precision)

	# Creating local variables :
	# --------------------------
	D = np.zeros(shape=len(N))   # plane constant
	M = np.zeros(shape=(3,3))    # system to solve

	# Creating output :
	# -----------------
	faces_ = []                     # the faces of the domain are stored here
	points_ = np.empty(shape=(0,3)) # all the points computed are stored here

	# Local parameters :
	# ------------------
	Esp  = 0.00001   # precision
	Esp2 = 0.00001	 
	doublon = False

	# Computing plane constants :
	# ---------------------------
	for k in range(0,len(D)):
		D[k] = sum(A[k,:]*N[k,:])

	# Loop for creating all intersections points :
	# --------------------------------------------
	for i in range(0,len(N)):
		for j in range(i+1,len(N)):
			for k in range(j+1,len(N)):
				M     = np.vstack((N[i,:],N[j,:],N[k,:]))
				delta = np.linalg.det(M)
				if abs(delta)>Esp:
					X = np.linalg.solve(M,[D[i],D[j],D[k]])
					# Check if X is in domain : 
					# -------------------------
					if (FinDomain(X,A,N,O) 
							and (P_boundaries[0])<=X[0]<=P_boundaries[1]
							and Q_boundaries[0]<=X[1]<=Q_boundaries[1]
							and U_boundaries[0]-0.0001<=X[2]<=U_boundaries[1]
							):
						doublon = False
						for l in range(0,len(points_)):
							if np.linalg.norm(X-points_[l,:])<Esp2:
								doublon = True
								break
						if(doublon!=True):
							points_ = np.vstack((points_,X))

	# Loop for creating domain faces :
	# --------------------------------
	for i in range(0,len(N)):
		f = Face(N[i,:])
		for j in range(0,len(points_)):
			if abs(sum(points_[j,:]*N[i,:])-D[i])<Esp:
				# Add point to face
				f.addPoint(points_[j,:])

		if len(f.points)>=3:
			# Order points of face
			f.orderPoints()
			# Add face to faces_
			faces_.append(f)

	return faces_, points_

def FdomainExport(faces, points, filename):
	"""
		Export domain in OFF format file
		Inputs :
			- faces    : composed of normal vector and ordered vertex
			- points   : contain all the vertex of the domain
			- filename : desired name for the output file
		Ouput  :
			- file in OFF format
		Used in :
			- 
	"""	
	# Computing number of edges using Euler's formula :
	# -------------------------------------------------
	edges = len(points) + len(faces) - 2

	# Writing ouput file :
	# --------------------
	with open(os.path.join(filename + ".off"),'w') as f:
		# Writing header :
		# -----------------
		f.write("OFF\n")
		f.write(repr(len(points)) + " " + repr(len(faces)) + " " + repr(edges) + "\n")
		# Writing points :
		# ----------------
		np.savetxt(f,points,fmt='%f')
		# Writing faces :
		# ---------------
		for i in range(0,len(faces)):
			# Writing number of points composing face i :
			# -------------------------------------------
			f.write(repr(faces[i].points.shape[0]) + " ")
			for j in xrange(faces[i].points.shape[0]):
				# Writing index of points composing face i :
				# ------------------------------------------
				np.savetxt(f,np.where(np.all(faces[i].points[j,:]==points,axis=1)),fmt='%i',newline=' ')
			f.write("\n")

def FconvexHullExport(points, filename):
	"""
		Create convex hull and export it to .off file
		Inputs :
			- points   : contain all the vertex of the domain
			- filename : desired name for the output file
		Ouput  :
			- file in OFF format
		Used in :
			- MachineStableDomainSearch.py
	"""	
	points_clean = np.zeros(shape=(0,3))
	convH = ConvexHull(points)

	edges = len(convH.vertices) + len(convH.simplices) - 2

	# Writing ouput file :
	# --------------------
	with open(os.path.join(filename + ".off"),'w') as f:
		# Writing header :
		# -----------------
		f.write("OFF\n")
		f.write(repr(len(convH.vertices)) + " " + repr(len(convH.simplices)) + " " + repr(edges) + "\n")
		# Writing points :
		# ----------------
		for vertice in convH.vertices:
			np.savetxt(f,convH.points[vertice,:],fmt='%f',newline=' ')
			f.write("\n")
			points_clean = np.vstack((points_clean,convH.points[vertice,:]))
		# Writing faces :
		# ---------------
		# np.savetxt(f,convH.simplices,fmt='%i')
		for i in range(0,len(convH.simplices)):
			# Writing number of points composing face i :
			# -------------------------------------------
			f.write(repr(len(convH.simplices[i])) + " ")
			for j in range(0,len(convH.simplices[i])):
				np.savetxt(f,np.where(np.all(convH.points[convH.simplices[i,j],:]==points_clean,axis=1)),fmt='%i',newline=' ')
			f.write("\n")

	compteur = 0
	doublon = False
	for i in range (0,len(convH.equations)):
		doublon = False
		for j in range(i+1,len(convH.equations)):
			if np.all(convH.equations[i] == convH.equations[j]):
				doublon = True
				break
		if doublon == False:
			compteur += 1

	# print "Nombre d'equations : ", len(convH.equations)
	# print "Nombre reel        : ", compteur

def Ffitfunction(X,a,b,c,d):
	"""
		function to fit to data
		Inputs :
			- X : (P,Q) Active power and reactive power
		Ouput  :
			- Z : U Tension
		Not used
	"""		
	x,y = X
	return a*x**4 + b*y**2 + c*x*y + d

def FexportCSV(points,filename,Machine_ref_name,Machine_int_name,U_nominal,O):
	"""
		Create an export of planes equation to be used by projector. There should not be similar equation
		twice. Comments are specified by # line started.
		- Inputs :
			- points   : the planes are computed from the convex hull of these points
			- filename
			- Machine_ref_name : name of the power generator given by the dictionnary
			- Machine_int_name : name of the power generator given in eurostag
			- U_nomilan : nominal tension of generator
			- O : stable origin 
		- Outputs : 
			- filename.txt
		- Used in : MachineStableDomaineSearch.py/main

	"""

	# Compute convexe hull of input points : 
	# --------------------------------------
	convH = ConvexHull(points)

	equation_reel = np.zeros(shape=(0,4))

	# Checking for duplicate :
	# ------------------------
	compteur = 0
	doublon = False
	for i in range (0,len(convH.equations)):
		doublon = False
		for j in range(i+1,len(convH.equations)):
			if np.all(convH.equations[i] == convH.equations[j]):
				doublon = True
				break
		if doublon == False:
			equation_reel = np.vstack((equation_reel,convH.equations[i]))
			compteur += 1

	# Checking if stable origin is inside domain :
	# --------------------------------------------
	for i in range(0,equation_reel.shape[0]):
		if (np.dot(O,equation_reel[i,:3]) > 0.00001 - equation_reel[i,3]):
			print "Problème avec le point d'origine"
			print "Equation :", equation_reel[i,:] 
			print "\n"
	# Checking if points of convex hull are inside domain :
	# -----------------------------------------------------
	for point in convH.points:
		for i in range(0,equation_reel.shape[0]):
			if (np.dot(point,equation_reel[i,:3]) > 0.00000001 - equation_reel[i,3] ):
				print point
				print equation_reel[i,:]
				print "Warning un point n'est pas dans le domaine !!"

	# Changing the last coefficient to have Ax < b
	equation_reel[:,3] = -equation_reel[:,3]

	# Writing ouput file :
	# --------------------
	with open(os.path.join(filename + ".txt"),'w') as f:

		# Writing domain info :
		# ---------------------
		f.write("# " + time.strftime("%d/%m/%y") + "\n")
		f.write("# P_max(MW)   : " + str(points[np.where(points[:,0] == np.max(points[:,0]))[0][0],:]) + "\n")
		f.write("# P_min(MW)   : " + str(points[np.where(points[:,0] == np.min(points[:,0]))[0][0],:]) + "\n")
		f.write("# Q_max(MVar) : " + str(points[np.where(points[:,1] == np.max(points[:,1]))[0][0],:]) + "\n")
		f.write("# Q_min(MVar) : " + str(points[np.where(points[:,1] == np.min(points[:,1]))[0][0],:]) + "\n")
		f.write("# U_max(kV)   : " + str(points[np.where(points[:,2] == np.max(points[:,2]))[0][0],:]) + "\n")
		f.write("# U_min(kV)   : " + str(points[np.where(points[:,2] == np.min(points[:,2]))[0][0],:]) + "\n")

		# Writing header :
		# ----------------
		f.write("#num id P(MW) Q(MVar) V(p.u) RHS(lt) " + "Vnominal(kV)" + " " + "id_internal" + "\n")
		# Writing planes :
		# ---------------
		for i in range(0,equation_reel.shape[0]):
			f.write(str(i+1) + " " + Machine_ref_name + " ")
			np.savetxt(f,equation_reel[i,:],fmt='%f',newline=" ")
			f.write(str(U_nominal) + " " + Machine_int_name + "\n")

def main(P_draw):
	"""
		Used to test the differents functions in this file
	"""
	p = 8
	A = np.zeros(shape=(p,3)) #store boundary points
	N = np.zeros(shape=(p,3)) #store boundary normal vector
	O = np.zeros(shape=(1,3))
	points = np.empty(shape=(0,3))

	P_boundaries = [0, 1]
	Q_boundaries = [-1, 1]
	U_boundaries = [0.7, 1.3]

	O[:] = [0.8,0,1]
	A[0,:] = [P_boundaries[0], rand.uniform(Q_boundaries[0],Q_boundaries[1]), rand.uniform(U_boundaries[0],U_boundaries[1])]
	A[1,:] = [P_boundaries[1], rand.uniform(Q_boundaries[0],Q_boundaries[1]), rand.uniform(U_boundaries[0],U_boundaries[1])]
	A[2,:] = [rand.uniform(P_boundaries[0],P_boundaries[1]), Q_boundaries[0], rand.uniform(U_boundaries[0],U_boundaries[1])]
	A[3,:] = [rand.uniform(P_boundaries[0],P_boundaries[1]), Q_boundaries[1], rand.uniform(U_boundaries[0],U_boundaries[1])]
	A[4,:] = [rand.uniform(P_boundaries[0],P_boundaries[1]), rand.uniform(Q_boundaries[0],Q_boundaries[1]), U_boundaries[0]]
	A[5,:] = [rand.uniform(P_boundaries[0],P_boundaries[1]), rand.uniform(Q_boundaries[0],Q_boundaries[1]), U_boundaries[1]]
	A[6,:] = [0.9*P_boundaries[1],0.9*Q_boundaries[1],0.9*U_boundaries[1]]
	A[7,:] = [P_boundaries[0]+0.1,0.8*Q_boundaries[0],1.1*U_boundaries[0]]
	
	N[0,:] = [-1, 0, 0]
	N[1,:] = [1, 0, 0]
	N[2,:] = [0, -1, 0]
	N[3,:] = [0, 1, 0]
	N[4,:] = [0, 0, -1]
	N[5,:] = [0, 0, 1]
	N[6,:] = [0.5,0.5,0.5]
	N[7,:] = [0.5,0.5,0.5]

	faces, points = FdomainCreation(A,N,P_boundaries,Q_boundaries,U_boundaries,O)
	# print "Point d'intersection : ", points 
	FdomainExport(faces, points, "test")
	FconvexHullExport(A,"convHtest")

	param, opt = curve_fit(Ffitfunction,A[:,:2].T,A[:,2])

	XX,YY = np.meshgrid(np.arange(0,1,0.01),np.arange(-1,1,0.02))
	print XX.shape[0], XX.shape[1]
	print YY.shape[0], YY.shape[1]
	ZZ = np.zeros(shape=(XX.shape[0],XX.shape[1]))

	for i in range(0,XX.shape[1]):
		PP = np.vstack((XX[:,i],YY[:,i]))
		print PP.shape[0], PP.shape[1]
		ZZ[:,i] = Ffitfunction(PP,param[0],param[1],param[2],param[3])

	# fig = plt.figure()
	# ax = fig.add_subplot(111, projection='3d')
	# ax.plot_surface(XX,YY,ZZ,rstride=5,cstride=5,alpha=0.5)
	# ax.scatter(A[:,0],A[:,1],A[:,2], c='r', s=50)
	# plt.xlabel('P')
	# plt.ylabel('Q')
	# ax.set_zlabel('U')
	# ax.axis('equal')
	# ax.axis('tight')
	# plt.show()
	FexportCSV(A,"testexport","TEST_True", "TEST_int", 10,O)
	
	print A

if __name__ == '__main__':

	P_draw = float(sys.argv[1])

	main(P_draw)	

















