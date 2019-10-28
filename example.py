#!/usr/bin/env python
# ~*~ coding:utf-8 ~*~

"""
Example file for module_helix_angle function
"""

import module_helix_angle as mha
import sys

################################################################################
#############Uncomment the command line that you want to test###################
################################################################################


# the limits list contain on each tuples, the first and the last atom of an 
# helix. If you want to compute angles for all helix in the protein, remove 
# "limits" in the command line.

limits = [(656,891),(1884,2194),(2382,2598)]


#-----------------Computes angles for a PDB file :-------------------------------

##Computes helix-helixBundle angles for a pdb file :
mha.PDBHelix_vs_helix("data/1H6I.pdb")

##Computes helix-helix angles for a pdb file :
# mha.PDBHelix_vs_helixBunble("data/1H6I.pdb")

##Computes helix-helix torsion angles for a pdb file :
# mha.PDBhelix_torsion_angle("data/1H6I.pdb")



#---------------Computes angles during a MD simulation:-------------------------

## CAUTION : if you to use the next functions you must have 1HI.gro
## and 1H6I_traj.xtc in data. You can dowload these files on this link:
## http://www.dsimb.inserm.fr/~santuz/projet_court_helices/

## For all the calculation during the MDsimulation, the results is write on an
## output file. If you don't want to write this file remove the "w" in the command
## line else put any character.

## Compute Helix-Helix angles during a MD simulation :
# results,name = mha.MDhelix_vs_helix("data/1H6I.gro","data/1H6I_traj.xtc","w",limits)
# mha.plot_all(results,name,"helix-helix")
# mha.plot_angle("H1-H2",results,test)

## Compute Helix-HelixBundle angles during a MD simulation :
## in this command line, the output_file is not write.
# results,name = mha.MDhelix_vs_helixbundle("data/1H6I.gro","data/1H6I_traj.xtc",0,limits)
# mha.plot_all(results,name,"helix-bundlehelix")

## Compute Helix-Helix torsion angles during a MD simulation :
## This calculation can take the limits argument, this fuction compute the torsion
#for helix who are in contact.
# results,name = mha.MDhelix_torsion_angle("data/1H6I.gro","data/1H6I_traj.xtc","w")
# mha.plot_all(results,name,"helix-helix_torsion_angle")