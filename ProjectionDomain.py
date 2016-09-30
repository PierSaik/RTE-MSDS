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
import numpy as np
import time
from scipy.spatial import ConvexHull
import multiprocessing as mp

import matplotlib.pyplot as plt
from matplotlib import collections as mc
# Script :
# --------

def FinDomain(point, domain_data):
	"""
		Check if the point is inside the domain
		Inputs : 
			- point
			- domain_raw 
		Outputs :
			- True/False
		Used in :
			- FcalcDistance

	"""	
	# Local Variable :
	# --------------- 
	eps = 0.0001   # Tolerance

	# Testing point : 
	# ---------------
	for i in range(0,domain_data.shape[0]):
		if (np.dot(domain_data[i,:3],point) - domain_data[i,3]>= eps):
			return False
			break
	
	return True

def FinPlane(point, planes):
	"""
		Check if the point is inside the domain
		Inputs : 
			- point
			- PQ_planes
		Outputs :
			- True/False
		Used in :
	"""	
	# Local Variable :
	# ----------------
	eps = 0.001    # precision

	# Testing point : 
	# ---------------
	for i in range(0,len(planes)):
		if (np.dot(planes[i,:2],point) + planes[i,2] >= eps):
			return False
			break
	
	return True

def FcalcDistance(gen_data, domain_data):
	"""
		Compute the distance from the input point to the domain.
		Inputs : 
			- pgen_data
			- domain_data 
		Outputs :
			- -1 -> if the point is inside the domain
			- dom_dist -> if the point is outsite the domain
		Used in :
			- main
	"""	

	dom_dist_PQ = 0.
	dom_dist_UQ = 0.
	point_proj_PQ_OK = np.zeros(shape=(2))
	point_proj_UQ_OK = np.zeros(shape=(2))

	if FinDomain(gen_data,domain_data):
		return -1., 0, 0, -1, 0, 0 
	else:
		points_corner_PQ, PQ_planes = FcomputePQdiag(domain_data, gen_data[2])

		for i in range(0,PQ_planes.shape[0]):
			d_plane = (np.dot(gen_data[:2],PQ_planes[i,:2]) + PQ_planes[i,2])/np.linalg.norm(PQ_planes[i,:2])
			if (d_plane > 0):
				point_projection = gen_data[:2] - d_plane*PQ_planes[i,:2] 
				if (FinPlane(point_projection, PQ_planes) and d_plane >= dom_dist_PQ):
					dom_dist_PQ = d_plane
					point_proj_PQ_OK = gen_data[:2] - d_plane*PQ_planes[i,:2] 
		for i in range(0,len(points_corner_PQ)):
			d_corner = np.linalg.norm(gen_data[:2] - points_corner_PQ[i,:])
			if (d_corner < dom_dist_PQ):
				dom_dist_PQ = d_corner
				point_proj_PQ_OK = points_corner_PQ[i,:2]

		points_corner_UQ, UQ_planes, status = FcomputeUQdiag(domain_data, gen_data[0])

		if (status != 0):
			for i in range(0,UQ_planes.shape[0]):
				d_plane = (np.dot(gen_data[1:],UQ_planes[i,:2]) + UQ_planes[i,2])/np.linalg.norm(UQ_planes[i,:2])
				if (d_plane > 0):
					point_projection = gen_data[1:] - d_plane*UQ_planes[i,:2] 
					if (FinPlane(point_projection, UQ_planes) and d_plane >= dom_dist_UQ):
						dom_dist_UQ = d_plane
						point_proj_UQ_OK = gen_data[1:] - d_plane*UQ_planes[i,:2] 
			for i in range(0,len(points_corner_UQ)):
				d_corner = np.linalg.norm(gen_data[1:] - points_corner_UQ[i,:])
				if (d_corner < dom_dist_UQ):
					dom_dist_UQ = d_corner
					point_proj_UQ_OK = points_corner_UQ[i,:2]
		else:
			dom_dist_UQ      = -99999
			points_corner_UQ = 0.
			point_proj_UQ_OK = 0.

		return dom_dist_PQ, points_corner_PQ, point_proj_PQ_OK, dom_dist_UQ, points_corner_UQ, point_proj_UQ_OK

def FcomputePQdiag(data_raw, U_draw):
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
			M_sys = np.vstack((data_raw[i,0:2],data_raw[j,0:2]))
			delta = np.linalg.det(M_sys)
			if (abs(delta) > eps):
				X = np.linalg.solve(M_sys,[data_raw[i,3] - data_raw[i,2]*U_draw, data_raw[j,3] - data_raw[j,2]*U_draw])
				# Check if X is in domain : 
				# -------------------------
				testeur = X
				testeur= np.append(X,U_draw)
				if (FinDomain(testeur, data_raw)):
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
	PQ_planes     = Hull.equations

	return points_sorted, PQ_planes

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
				testeur = X
				testeur= np.append(P_draw,X)
				if (FinDomain(testeur, data_raw)):
					doublon = False
					for l in range(0,len(points)):
						if np.linalg.norm(X-points[l,:])<eps2:
							doublon = True
							break
					if(doublon!=True):
						points = np.vstack((points,X))

	# Computing Hull of points to order them :
	# ----------------------------------------
	if points.shape[0] != 0:
		Hull = ConvexHull(points)
		points_sorted = points[Hull.vertices,:]
		UQ_planes     = Hull.equations
		status = 1
	else:
		points_sorted = 0
		UQ_planes = 0
		status = 0


	return points_sorted, UQ_planes, status

def FfindMatch(gen_id, ampl_domain_id, ampl_domain_data, ampl_domain_Unom):
	"""
		Find the matching domain data of the generator in input..
		Inputs : 
			- gen_id : id of generator
			- ampl_domain_id : list of all the occurence of generators id
			- ampl_domain_data : list of all the planes of generators
			- ampl_domain_Unom : list of all the Unom of generators
		Outputs :
			- domain_data : planes of the gen_id
			- Unom : nominal tension of gen_id
			- 0 : if generator is not found
		Used in :
			- main
	"""	

	domain_index = ampl_domain_id==gen_id
	if any(domain_index==True) :
		# print "domain index : ", domain_index
		domain_data = ampl_domain_data[domain_index==True,:]
		Unom = ampl_domain_Unom[domain_index==True][0]

		return domain_data, Unom

	else:
		return 0, 0
	
def FProjectorCheck(proj_folder):
	"""
		Compute the distance between the ordered point of generator and the domain for each 
		machine in the situation. 
		Inputs : 
			- ampl_net_gen_path : path of the ordered points
			- ampl_domain_path  : path of the generators domain
		Outputs :
			- network_to_domain.txt : File containing the distance for each generators
		Used in :
			- main
	"""
	time_start = time.clock()
	fig, ax = plt.subplots()

	gen_info_path = os.path.join('generators_in_info.txt')
	gen_info_data = np.loadtxt(gen_info_path, comments='#', dtype=str)

	print "Currently Working on : ", proj_folder, '...'
	ampl_net_gen_path = os.path.join(proj_folder, "ampl_network_generators.txt")
	ampl_domain_path  = os.path.join(proj_folder, "ampl_generators_domains.txt")

	# Results :
	# ---------
	net_dom_dist_PQ = np.zeros(shape=(0,2))
	net_dom_dist_UQ = np.zeros(shape=(0,2))

	# Extracting data :
	# -----------------
	ampl_net_gen_data_raw = np.loadtxt(ampl_net_gen_path, comments='#', usecols=[11,12,13])
	ampl_net_gen_id_raw   = np.loadtxt(ampl_net_gen_path, comments='#', usecols=[16], dtype=str)
	ampl_domain_data      = np.loadtxt(ampl_domain_path, comments='#', usecols=[2,3,4,5] )
	ampl_domain_id        = np.loadtxt(ampl_domain_path, comments='#', usecols=[1] ,dtype=str)
	ampl_domain_Unom      = np.loadtxt(ampl_domain_path, comments='#', usecols=[6])
	# projector_results     = np.loadtxt(projector_results_path, comments='#' , usecols=[1,2,3])
	# projector_id          = np.loadtxt(projector_results_path, comments='#', usecols=[0], dtype=int)

	# Cleaning network data :
	# -----------------------
	for i in range(0,ampl_net_gen_id_raw.shape[0]):
		ampl_net_gen_id_raw[i] = ampl_net_gen_id_raw[i].strip('"')
		ampl_net_gen_data_raw[i,:] = [ampl_net_gen_data_raw[i,1],ampl_net_gen_data_raw[i,2],ampl_net_gen_data_raw[i,0]]

	# Computing distance for each generator in list :
	# -----------------------------------------------
	for i, gen_id in enumerate(ampl_net_gen_id_raw):
		if (ampl_net_gen_data_raw[i,0]!=-1 and ampl_net_gen_data_raw[i,1]>0.):

			gen_data = ampl_net_gen_data_raw[i,:]
			domain_data, Unom = FfindMatch(gen_id, ampl_domain_id, ampl_domain_data, ampl_domain_Unom)

			if (Unom != 0):

				Snom = float(gen_info_data[gen_id == gen_info_data[:,0],2][0])
				gen_data[2] *= Unom

				dom_dist_PQ, points_corner_PQ, points_proj_PQ_OK, dom_dist_UQ, points_corner_UQ, points_proj_UQ_OK = FcalcDistance(gen_data, domain_data)
				if dom_dist_PQ >0.:
					dom_dist_PQ = dom_dist_PQ/Snom
					
				net_dom_dist_PQ = np.vstack((net_dom_dist_PQ,[gen_id, dom_dist_PQ]))
				net_dom_dist_UQ = np.vstack((net_dom_dist_UQ,[gen_id, dom_dist_UQ]))

				if dom_dist_PQ > 0.:

					# gen_proj = projector_results[projector_id==(i+1),:]
					# gen_proj = [gen_proj[0,1], gen_proj[0,2] , gen_proj[0,0]]
					plt.fill(points_corner_PQ[:,1], points_corner_PQ[:,0], color='purple', alpha=0.5, label='domain')
					plt.scatter(gen_data[1],gen_data[0], color='red' , s=50 , label='set point')
					plt.scatter(points_proj_PQ_OK[1], points_proj_PQ_OK[0], color='green' , s=50, label='proj')
					lc = mc.LineCollection([[(gen_data[1],gen_data[0]),(points_proj_PQ_OK[1],points_proj_PQ_OK[0])]], colors='red', linewidths=2)
					# plt.scatter(gen_proj[0],gen_proj[1], color='green' , s=50 , label='projection')
					ax.add_collection(lc)
					plt.ylabel('P Active Power (MW)', fontsize='x-large')
					plt.xlabel('Q Reactive Power (MVar)', fontsize='x-large')
					plot_title = "Diagramme PQ a " + str(round(gen_data[2],1)) + " kV pour " + gen_id
					plt.title(plot_title, fontsize='xx-large')
					plt.legend(loc='lower left', shadow=True, fontsize='xx-large', scatterpoints = 1)
					fig_name = os.path.join(proj_folder, gen_id + '_PQ.png')
					fig.savefig(fig_name)
					fig.clear()
					ax.clear()

				if dom_dist_UQ > 0.:

					# gen_proj = projector_results[projector_id==(i+1),:]
					# gen_proj = [gen_proj[0,1], gen_proj[0,2] , gen_proj[0,0]]
					plt.fill(points_corner_UQ[:,1], points_corner_UQ[:,0], color='purple', alpha=0.5, label='domain')
					plt.scatter(gen_data[2],gen_data[1], color='red' , s=50 , label='set point')
					plt.scatter(points_proj_UQ_OK[1], points_proj_UQ_OK[0], color='green' , s=50, label='proj')
					lc = mc.LineCollection([[(gen_data[2],gen_data[1]),(points_proj_UQ_OK[1],points_proj_UQ_OK[0])]], colors='red', linewidths=2)
					# plt.scatter(gen_proj[0],gen_proj[1], color='green' , s=50 , label='projection')
					ax.add_collection(lc)
					plt.ylabel('Q Reactive Power (MVar)', fontsize='x-large')
					plt.xlabel('U Tension (kV)', fontsize='x-large')
					plot_title = "Diagramme UQ a " + str(round(gen_data[0],1)) + " MW pour " + gen_id
					plt.title(plot_title, fontsize='xx-large')
					plt.legend(loc='lower left', shadow=True, fontsize='xx-large', scatterpoints = 1)
					fig_name = os.path.join(proj_folder, gen_id + '_UQ.png')
					fig.savefig(fig_name)
					fig.clear()
					ax.clear()

			else:
				# No Data
				net_dom_dist_PQ = np.vstack((net_dom_dist_PQ,[gen_id,"-99999"]))
				net_dom_dist_UQ = np.vstack((net_dom_dist_UQ,[gen_id,"-99999"]))
		else:
			# Not connected
			net_dom_dist_PQ = np.vstack((net_dom_dist_PQ,[gen_id,"-100"]))
			net_dom_dist_UQ = np.vstack((net_dom_dist_UQ,[gen_id,"-100"]))

	time_elapsed = (time.clock() - time_start)
	print "Done with : ", proj_folder, "in : ", time_elapsed, "s"
	plt.close(fig)

	net_dom_dist_PQ = np.vstack((["gend_id", proj_folder],net_dom_dist_PQ))
	net_dom_dist_UQ = np.vstack((["gend_id", proj_folder],net_dom_dist_UQ))


	return net_dom_dist_PQ, net_dom_dist_UQ

def main():

	net_dom_resultsPQ = np.zeros(shape=(0,0))
	net_dom_resultsUQ = np.zeros(shape=(0,0))
	net_dom_results_rawPQ = []
	net_dom_results_rawUQ = []
	proj_folder_list = []

	# Get list of worker from cpu count
	mp.freeze_support()
	PROCESSES = mp.cpu_count()-1

	for proj_folder in os.listdir('.'):
		if proj_folder.startswith('itesla_projector_'):
			proj_folder_list.append(proj_folder)

	# A, B =  FProjectorCheck(proj_folder_list[0])
	
	pool = mp.Pool(processes=PROCESSES)
	net_dom_results_raw = list(pool.imap_unordered(FProjectorCheck, proj_folder_list))
	time.sleep(0.1)
	pool.close()
	pool.join()
	# print net_dom_results_raw


	for i in range(0,len(net_dom_results_raw)):
		net_dom_results_rawPQ.append(net_dom_results_raw[i][0])
		net_dom_results_rawUQ.append(net_dom_results_raw[i][1])

	# np.savetxt('net_dom_results_rawPQ.txt', net_dom_results_rawPQ, fmt="%s")
	# np.savetxt('net_dom_results_rawUQ.txt', net_dom_results_rawUQ, fmt="%s")	

	for j in range(0,len(net_dom_results_rawPQ)):
		if (net_dom_resultsPQ.shape[0] == 0):
			net_dom_resultsPQ = net_dom_results_rawPQ[j]
		else:
			net_dom_resultsPQ = np.hstack((net_dom_resultsPQ, np.zeros((net_dom_results_rawPQ[j].shape[0],1))))
			for i in range(0,net_dom_resultsPQ.shape[0]):	
				net_dom_resultsPQ[i,-1] = net_dom_results_rawPQ[j][net_dom_resultsPQ[i,0]==net_dom_results_rawPQ[j][:,0],1][0]

	for j in range(0,len(net_dom_results_rawUQ)):
		if (net_dom_resultsUQ.shape[0] == 0):
			net_dom_resultsUQ = net_dom_results_rawUQ[j]
		else:
			net_dom_resultsUQ = np.hstack((net_dom_resultsUQ, np.zeros((net_dom_results_rawUQ[j].shape[0],1))))
			for i in range(0,net_dom_resultsUQ.shape[0]):	
				net_dom_resultsUQ[i,-1] = net_dom_results_rawUQ[j][net_dom_resultsUQ[i,0]==net_dom_results_rawUQ[j][:,0],1][0]

	np.savetxt('domain_to_network_dist_PQ.txt', net_dom_resultsPQ, fmt="%s")
	np.savetxt('domain_to_network_dist_UQ.txt', net_dom_resultsUQ, fmt="%s")

if __name__=='__main__':

	# ampl_net_gen_path = sys.argv[1]     # version of progrom used
	# ampl_domain_path  = sys.argv[2]     # visualisation on/off
	# projector_results_path = sys.argv[3]

	main()
