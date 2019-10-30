# HelixAngle
 Python module: Calculation of helices angle during a MD

![](pictures/anglelogo.png "logo" )

 ----------------------------------------------------------------------------
    		module_helix_angle.py (16-Sept-2015) 	version 1.00   
 ----------------------------------------------------------------------------
 
 
 Copyright (c) Emmanuel Edouard MOUTOUSSAMY M2BI

 Contents
----------
 1. Introduction
 2. Necessary packages and software
 3. Examples
 4. Contact


#1. Introduction
-------------------
This packages allow to determine differents angles between alpha-helix:
	- Helix - helix angle
	- Helix - Helix bundle angle
	- helix - Helix torsion angle

The helix position can be given by the user or calculated with the dssp program.
The computing time depends of the size of the protein.

This packages contain 6 functions :

1. PDBHelix_vs_helix(top,traj,output = 0,limits = 0):
	Compute helix -helix angle for a pdb file. A topology and a trajectory file
	are obligatory. If you want to write the results in an outpufile, the output
	argument must be "w". And if you want to give the helix position, limits must
	be a list which contain the helix position (ex. limits = [(12,45),(56,75)]).

2. PDBHelix_vs_helixBundle(top,traj,output = 0,limits = 0):
	Compute helix - helixBundle angle for a pdb file. The arguments is the same
	than 1.

3. PDBHelix_torsion_angle(top,traj,output = 0):
	Compute helix - helixBundle angle for a pdb file. 

4. MDHelix_vs_helix(top,traj,output = 0,limits = 0):
	Compute Helix-helix angle during a MD simulation. If output = "w", the
	results is write in an output file.

5. MDHelix_vs_helixBundle(top,traj,output = 0,limits = 0):
	Compute Helix-helixBundle angle during a MD simulation. The arguments is the
	same than 4.

6. MDHelix_torsion_angle(top,traj,output = 0):
	Compute Helix torsion angle during a MD simulation. The arguments is the same
	than 4.

	
 #2. Necessary packages and software
-----------------------------------

The necessary packages are :
modules python:
-SYS (basic python package)
-OS (basic python package)

-NUMPY (http://www.numpy.org):
Debian / Ubuntu : sudo apt-get install python-numpy
Mac OSX-MAVERICK (via MacPorts) : sudo port install py27-numpy

-MATPLOTLIB (http://matplotlib.org):
Debian / Ubuntu : sudo apt-get install python-matplotlib
Mac OSX-MAVERICK (via MacPorts) : sudo port install py27-matplotlib

-MDAnalysis :
Debian / Ubuntu : pip install MDAnalysis
Mac OSX-MAVERICK : pip install MDAnalysis

!Caution : the mkdssp program is essential for some options !!!!
http://pkgs.org/download/mkDSSP

 #3. Example
-----------------------------------
the example file contain differents uses of this program. For run example.py, 
Uncomment the lines that your want to test and type :

./example.py

The data file contain :
- a pdb file : 1H6I.pdb


 #4. Contact
---------------------------------
For all questions or informations:
e.e.moutoussamy@gmail.com
--------------------------------------------------------------------------------
