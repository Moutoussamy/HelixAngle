#!/usr/bin/env python
# ~*~ coding:utf-8 ~*~


"""
This script permit to calculate differetns angles between alpha-helix from a
PDB file or a MD simulation (CHARMM or GROMACS).
The differents angles are :
	- helix - helix angle
	- helix - helix bundle angle
	- helix - helix torsion angle
"""

import lib.function_helix_angle as fct
import lib.math_angle as ma
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
import os

__author__ = "Emmanuel Edouard MOUTOUSSAMY"
__version__  = "1.0.0"
__date__ = "2015/09"
__copyright__ = "CC_by_SA"
__dependencies__ = "os, sys, numpy, math, matplotlib, MDAnalysis"


#-----------------------------PDB-----------------------------------------------

def PDBHelix_vs_helix(pdb,limits_file = 0):
	"""
	Compute the angle between for all combination of two diffrents helix.
	all the results is write on an output file.
	Args : a pdb file and a limits file (optional)
	"""
	if limits_file == 0:
		ss = fct.secondary_structure(pdb) # compute the helix position
	else:
		ss = fct.get_ss(limits_file,pdb)# get the helix position from helix list

	list_helix = fct.search_helix(pdb,ss) # extract helix from a pdb
	fct.print_limits(list_helix) # print and write the helix position

	fct.get_inertie_mtx(list_helix) # compute the inertia matrix
	fct.helix_evec(list_helix) # compute the eigen vector for each helix

	print "\n"
	print "Results :"
	output = open("helix_helix_angle.dat","w")
	header = "#Helix angle(degrees)\n"
	output.write(header) # write the header on the output file
	print header

	start = 1
	# this double loop allows to compute angle between each helix couple
	for i in range(len(list_helix)):
		axis1 = list_helix[i].eig_vec[:,0] # axis of helix 1
		for j in range(start,len(list_helix)):
			axis2 = list_helix[j].eig_vec[:,0] # axis for helix 2
			axis1 = fct.check_vec_sens(axis1,axis2) # check the vectors direction
			angle = math.acos(axis1.dot(axis2))*57.3
			out = "%s-%s %.3f"%(list_helix[i].name,list_helix[j].name,angle)
			output.write("%s\n"%out) # write the angle on the output file
			print out
		start += 1
	print "The resutls file : helix_helix_angle.dat"
	print ""
	output.close()

def PDBHelix_vs_helixBunble(pdb,limits_file = 0):
	"""
	Allows to compute the angle between an helix and the protein axis on a pdb
	file.
	Args : a pdb file and a limits file (optional)
	Return : a file which contain all the angle calculation
	"""
	if limits_file == 0:
		ss = fct.secondary_structure(pdb) # compute the helix position
	else:
		ss = fct.get_ss(limits_file,pdb)# get the helix position from helix list

	list_helix = fct.search_helix(pdb,ss) # extract helix from a pdb
	fct.print_limits(list_helix) # print and write the helix position

	fct.get_inertie_mtx(list_helix) # compute the inertia matrix
	fct.helix_evec(list_helix)# compute the eigen vector for each helix
	hb_axis = fct.helix_bundle_evec(list_helix) # compute the helix bundle axis
	output = open("prot_angle.dat","w")
	output.write("#Helix angle(degrees)\n")
	axis1 = hb_axis[:,0]
	i = 1
	print "\nResults : "
	for helix in list_helix:
		axis2 = helix.eig_vec[:,0]
		#Compute the helix-helix bundle angle :
		angle = math.acos(axis1.dot(axis2))*57.3
		out = "%s-helix_bundle_axis %.3f\n"%(helix.name,angle)
		output.write(out)
		print out
		i += 1
	output.close()

def PDBhelix_torsion_angle(pdb,limits_file = 0):
	"""
	Compute the torsion angle for helix which is contact in a pdbfile
	"""
	if limits_file == 0:
		ss = fct.secondary_structure(pdb) # compute the helix position
	else:
		ss = fct.get_ss(limits_file,pdb)# get the helix position from helix list

	list_helix = fct.search_helix(pdb,ss) # extract helix from a pdb
	fct.print_limits(list_helix) # print and write the helix position
	contact = fct.WhichIsInContact(list_helix)# Compute which helix is in contact
	fct.get_inertie_mtx(list_helix) # compute the inertia matrix
	fct.helix_evec(list_helix) # Compute the eigen vectors for each helix
	output = open("torsion_angle.dat","w")
	header = "#Helix torsion_angle(degrees)\n"
	output.write(header)
	print header
	start = 1
	for i in range(len(list_helix)):
		h1 = list_helix[i]
		for j in range(start,len(list_helix)):
			h2 = list_helix[j]
			ang = "%s-%s"%(h1.name,h2.name)
			if ang in contact:
				# Compute the torsion angle :
				angle = ma.compute_torsion_angle(h1,h2)
				out = "%s-%s %.3f\n"%(h1.name,h2.name,angle)
				output.write(out)
				print out
		start += 1
	output.close()

#-----------------------------MD simulaion--------------------------------------


def MDhelix_vs_helix(top,traj,output = 0,limits = 0):
	"""
	Compute the helix - helix angle during a MD simulation.
	If output = "w", an output file is write.
	Args : a topology, a trajectory file and an output (optional)
	Return : a numpy matrix which contain all angle and a list which contain all
	angle couple.
	"""
	#mk a directory and compute the number of frame in the trajectory:
	last_frame,dir_name = fct.get_pdb(top,traj) #number of frame

	row = 0
	if output == "w":
		output = open('results_helix_vs_helix.dat','w')

	for i in range(last_frame):
		print 'Compute helix-helix angle for frame : %i'%(i+1)
		if i == 0:
			if limits == 0:
				# Compute the helix positions:
				ss = fct.secondary_structure("%s/dssp_input.pdb"%dir_name)
			else :
				# Extract the helix positions from an input file
				ss = fct.get_ss(limits,"%s/dssp_input.pdb"%dir_name)

			#Extract helix from the frame:
			list_helix = fct.search_helix("%s/frame_1.pdb"%dir_name,ss)

			if i == 0:
				#Print and write the helix positions
				fct.print_limits(list_helix,"w")

			# Create a numpy matrix with the good dimensions:
			results_mtx = fct.mk_results_mtx(last_frame,list_helix)
			names = fct.header_hvsh(output,list_helix)

		# Extract helix from the frame:
		list_helix = fct.search_helix("%s/frame_%i.pdb"%(dir_name,i+1),ss)
		# Compute the inertia matrix:
		fct.get_inertie_mtx(list_helix)
		#Compute the eigen vectors :
		fct.helix_evec(list_helix)

		results_mtx,row = fct.compute_hvsh(i,output,list_helix,results_mtx,row)

	if output != 0:
		output.close()

	# os.system("rm -rf %s"%dir_name)
	return results_mtx,names

def MDhelix_vs_helixbundle(top,traj,output = 0,limits = 0):
	"""
	Compute the helix - helix bundle angle during a MD simulation.
	If output = "w", an output file is write.
	Args : a topology, a trajectory file and an output (optional)
	Return : a numpy matrix which contain all angle and a list which contain all
	the combination of angle.
	"""
	#mk a directory and compute the number of frame in the trajectory:
	last_frame,dir_name = fct.get_pdb(top,traj)
	row = 0
	if output == "w":
		output = open('results_helix_vs_helixbundle.dat','w')
	for i in range(last_frame):
		print 'Compute helix - helix bundle angle for frame : %i'%(i+1)
		if i == 0:
			if limits == 0 :
				# Compute the helix positions:
				ss = fct.secondary_structure("%s/dssp_input.pdb"%dir_name)

			else:
				# Extract the helix positions from an input file:
				ss = fct.get_ss(limits,"%s/dssp_input.pdb"%dir_name)
			list_helix = fct.search_helix("%s/frame_1.pdb"%dir_name,ss)
			for helix in list_helix:
				print helix.name,helix.res_start,helix.res_finish

			if i == 0:
				fct.print_limits(list_helix,"w")
			results_mtx = np.zeros((last_frame,len(list_helix)))
			names = fct.header_hvshb(output,list_helix)

		list_helix = fct.search_helix("%s/frame_%i.pdb"%(dir_name,i+1),ss)
		#Compute the inertia matrix :
		fct.get_inertie_mtx(list_helix)

		#Compute the eigen vectors :
		fct.helix_evec(list_helix)

		# Compute the helixbundle axis:
		baxe = fct.helix_bundle_evec(list_helix)
		results_mtx,row = fct.compute_hvshb(i,output,list_helix,baxe[:,0],\
		results_mtx,row)

	if output != 0:
		output.close()

	os.system("rm -rf %s"%dir_name)
	return results_mtx,names

def MDhelix_torsion_angle(top,traj,output = 0):
	#mk a directory and compute the number of frame in the trajectory:
	last_frame,dirname = fct.get_pdb(top,traj)
	row = 0
	if output == "w":
		output = open('results_helix_torsion_angle.dat','w')

	for i in range(last_frame):
		print 'Compute helix - helix torsion angle for frame : %i'%(i+1)

		if i == 0:
			ss = fct.secondary_structure("%s/dssp_input.pdb"%dirname)

		list_helix = fct.search_helix("%s/frame_%i.pdb"%(dirname,i+1),ss)

		if i == 0:
			fct.print_limits(list_helix,"w")
			contact = fct.WhichIsInContact(list_helix)
			results_mtx = np.zeros((last_frame,len(contact)))

			if output != 0:
				fct.header_hta(output,contact)

		fct.get_inertie_mtx(list_helix) #Compute inertia matrix
		fct.helix_evec(list_helix) # Compute eigen vectors
		fct.evec_correction(list_helix) # Correct the vector direction
		fct.get_mass_center(list_helix) # Compute the mass center
		results_mtx,row = fct.compute_hta(i,output,list_helix,\
		contact,results_mtx,row)

	if output != 0:
		output.close()

	os.system("rm -rf %s"%dirname)

	return results_mtx,contact

#----------------------------PLOT-----------------------------------------------

def plot_all(results_mtx,names,dir_name = "plot"):
	"""
	Write all the graphics.
	Args : the results matrix, a list which contain the name of the different
	angle and a directory name.
	"""
	dir_name = fct.mk_directory(dir_name) # Make a directory
	col = results_mtx.shape[1]
	for i in range(col):
		print "Writting graphic %s"%names[i]
		val = results_mtx[:,i]
		mean = np.mean(val)
		plt.plot(range(len(val)),val,'.',color = "k") # plot
		plt.title("%s (mean angle = %.2f degree)"%(names[i],mean)) # add title
		plt.ylabel("angle (degrees)") # add y label
		plt.xlabel("Time step") # add x label
		plt.axhline(y=mean, color = 'red') # add a line for the mean value
		plt.ylim((min(val)-10,max(val)+10))
		plt.savefig("%s/plot_%s.png"%(dir_name,names[i])) # save
		plt.clf()

def plot_angle(angle,names,results_mtx):
	"""
	Write just one angle.
	Args : an angle name such as 'H1-H2' and a list which contain the name of the
	different angle and the results matrix.
	"""
	col = "null"
	# This loop allow to find the col of the angle in the results matrix
	for i in range(len(names)):
		if names[i] == angle:
			col = i

	if col == "null":
		sys.exit("ERROR : angle not find")

	val = results_mtx[:,col]
	mean = np.mean(val)
	plt.plot(range(len(val)),val,'.',color = "k") #Plot
	plt.title("%s (mean angle = %.2f degree)"%(names[i],mean)) #Title
	plt.ylabel("angle (degrees)") # add y label
	plt.xlabel("Time step") # add x lab
	plt.axhline(y=mean, color = 'red') # add a line for the mean value
	plt.ylim((min(val)-10,max(val)+10))
	plt.savefig("plot_%s.png"%(angle))# save
	plt.clf()
