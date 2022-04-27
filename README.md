To the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche.

# goMD, a tool for the analysis of MD trajectories.


This program makes extensive use of the goChem Computational Chemistry library.
If you use this program, we kindly ask you to support it by to citing the library as:

R. Mera-Adasme, G. Savasci and J. Pesonen, "goChem, a library for computational chemistry", http://www.gochem.org.

The gonum library (http://www.gonum.org/) is also much used (and appreciated).


goMD is a program that can calculate  different parameters for each frame of a trajectory, printing the results to the standard output, to be plot against the frame number.
Plotting scripts are also included.

## Binary install

Download the appropiate binary and add its location to PATH. 
If you wish to use the xtc-enabled binary, you need to have the xdrfile 1.1 library from Gromacs installed (i.e. its location must be added to LD_LIBRARY_PATH).

##  Compilation
goMD works with the current release version of goChem. It also requires the silly chemical utils (github.com/rmera/scu) and, if XTC file format support is needed, the xdrfile library from Gromacs (www.gromacs.com). As it uses modules, you can simply use the following command:

```
go install github.com/rmera/gomd
```

To compile and install the program (compiling with XTC-trajectories support requires the xdrfile library installed, and the "xtc" compilation tag).

## Use

```
program [-skip=X -begin=Y] Task molfilename trajname task_specific_arguments
```

The format of both the moleculefile (PDB, GRO and XYZ are supported) as well as that of the trajectory (DCD, multiPDB, multiXYZ,old-AMBER,STS and XTC are supported, XTC requires the Xdrlibrary), are determined by the respective file's extension. 

### Tasks

goMD can perform several tasks, and it was designed so that the implementation of new tasks is rather simple. The task name is *not* case-sensitive. The current tasks and their flags are:


* Ramachandran:
```
	./gomd  [-skip=X -begin=Y ] Ramachandran pdbname xtcname "residuename1 chain1" "residuename2 chain2" ... "residuenameN chainN" [RGB total_frames]
```

	Will plot the phi and psi angles for each residue/chain pair given on each frame on the trajectory. If the RGB keyword and the number of total frames that will be read from the trajectory are given, in addition to the phi/psi pair, three more columns will be printed for each residue/chain pair: The RGB numbers for a color (0 to 255) which will progress from red to purple (the standard hue circle) with increasing frames. This is useful to plot each phi/psi with its own color and follow the proggression of them with the trajectory. A little Python tool for plotting these reults is included.

* RMSD: Plots the RMSD of the given selections against the coordinates in the reference PDB file for those selections.

```	
	./gomd [-skip=X -begin=Y] RMSD pdbname xtcname "selection1" "selection2" ... "selectionN"
```

* ClosestN: Given a selection, plots the distances for the closest N of an also given list of residue names to the selection, for each frame.

```
	./gomd [-skip=X -begin=Y] ClosestN pdbname xtcname "selection" "residuename1 residuename2 residuenameM" N
```

	Where N is a positive integer. Residue names are given in 3 letter format (ASP, HIS, etc). The name for the water molecules will vary with the force-field used (SOL, HOH, WAT, etc). Distances of a residue with itself will not be considered.   


* WithinCutoff: Similar to the previous but returns the number of residues of the type residuename1 or residuename2, etc within R Angstroms of selection in each frame.

```
	./gomd [-skip=X -begin=Y] WhithinCutoff pdbname xtcname "selection" "residuename1 residuename2 residuenameM" R
```

	Where R is a positive float. The rest of the syntax is equivalent to that for ClosestN.


* RDF: Similar to the previous but returns the RDF or MDDF for residues of the type residuename1 or residuename2, with respecto to the selection.

```
	./gomd [-skip=X -begin=Y] RDF  pdbname xtcname "selection" "residuename1 residuename2 residuenameM" 
```

* Distance: Plots the distance between pairs of selections for each frame of the trajectory. 

```
	./gomd [-skip=X -begin=Y] distance pdbname xtcname "selection1" "selection2" ... "selectionN" 
```
	
	Where N is an even number. The distances between selections 1 and 2, ... N-1 and N for each frame will be printed to stdout. If a selection has more than one atom, the center of mass for that selection is used.

* Angle: Plots the angle ABC between sets of 3 selections "A", "B" and "C", for each frame of the trajectory. 

```
	./gomd [-skip=X -begin=Y] angle pdbname xtcname "selection1" "selection2" ... "selectionN" 
```

	Where N divisible by 3.  If a selection has more than one atom, the center of mass for that selection is used.


* Shape: Plots the Planarity (oblate distortion) an Elongation (prolate distorion) indicators for the selections.

```
	./gomd [-skip=X -begin=Y ] Shape pdbname xtcname "selection1" "selection2" ... "selectionN" 
```
* PlanesAngle: For every 2 selections, calculate the best plane passing through the atoms of each selection, and then returns the angles (in degrees) between the normal vector to each plane.

```
	./gomd [-skip=X -begin=Y ] planesangle pdbname xtcname "selection1" "selection2" ... "selection(N-1)" "selectionN"
```
Where N is an even number.

* FixGMX: Will print a "fixed" version a Gromacs PDB, with the chains restored, and exit. The chain identifiers will be added in the following way: The first lines will be assigned chain "A". Whenever the residue number of a residue is smaller than the residue number of a previous residue, a new chain will be assumed, which will be assigned the chain "B", and so on with "C", "D", etc.

```
	./gomd fixgmx pdbname
```

* Super Superimposes the whole trajectory to the reference structure considering the atoms in selection to calculate the superposition. In order to employ the LOVO procedure, instead of defining a selection, use the *-lovo* option (see below).

```
	./gomd super pdbname xtcname "selection"
```



* Average: Prints the average structure over the trajectory
```
	./gomd [-skip=X -begin=Y ] average pdbname xtcname 
```

* InterByRes: Takes 2 selections and, for each frame, returns the ID for the residues in each selection that are part of the interface between both selections in that frame. The residue IDs are returned as floats, but they are, of course, integer numbers. If both selections belong to different chains, the residues for the second selections will be given as negative numbers, to avoid mixing up residues from different chains but the same residue ID.

```
	./gomd [-skip=X -begin=Y ] interByRes pdbname xtcname "selection1" "selection2"
```

* LOVO superposition calculation: While not strictly a task, goMD can employ LOVO ( 10.1371/journal.pone.0119264\n, 10.1186/1471-2105-8-306) to find the optimal atoms from a selection to superimpose. LOVO is invoked with the option 

```
./gomd -lovo N task pdbname xtcname "selection"
```

Where N is a number larger than 0, and will become the number of frames skipped during the LOVO procedure. The task will be performed afterwards. Reasonable task options are "stop" (does nothing after printing the LOVO results) and "super" (which will superimpose the whole trajectory based on the index determined in the LOVO procedure).
	
### Selections: 

The selections are defined in the following way: "RESID1,RESID2,RESID3-RESID3+N,RESIDN CHAIN ATNAME1,ATNAME2"
RESID are residue numbers. They can be separated by commas or, to specify a range, with dashes: 12,13,120,125-128,145  The former would select the residues 12,13,120,125,126,127,128,145
CHAIN must be a chain identifier such as "A". If chain is "ALL", every chain will be used.
ATNAME is a PDB atom name such as CA (alpha carbon). Hydrogen names may vary with the forcefield employed. if ALL is given, as the first atom name, all atoms in the selected residues will be consiered.


### Plotting the results

A few Python scripts are included to help with the result visualization.

* ramachandran.py is a simple script that will plot one set of Ramachandran angles (the ones in the columns 2 and 3 of the output file from gomd) along the trajectory, coloring the value for each frame in a different hue, so any conformationa change can be observed.

* interfaces.py plots the results of the InterByRes, as the fraction of the trajectory for which each residue in the selection is part of the interface.

* plots.py is a rather complete plotting tool for all the other plotting tasks. It is easy to use (for the simplest case of only 1 column of results, or 2 columns in total, it only requires the name of the file) but has a fair bit of options to customize the visualization.  It can plot only the selected columns, set the scale for the y-axis, plot an histogram for the results instead of the default results vs time plot, and more. Plots.py will ignore lines starting with "\#", "@" and "&", which allows it to plot Gromacs-produced xvg files.

* rdf2gomd.py is a small program to put the results of the rdf task in a regular by-column goMD format, which can then be plotted with plots.py. It supports putting the results of several rdf runs as columns in the same file.



goMD, a little program for the analysis of molecular dynamics simulations.
Copyright (c) 2017  Raul Mera Adasme.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.



The author thanks ANID for financial support under Proyecto Fondecyt N. 11160032 and N. 1200200
