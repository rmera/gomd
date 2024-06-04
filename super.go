/*

goMD, a little tool for the analysis of MD trajectories.


This program makes extensive use of the goChem Computational Chemistry library.
If you use this program, we kindly ask you support it by to citing the library as:

R. Mera-Adasme, G. Savasci and J. Pesonen, "goChem, a library for computational chemistry", http://www.gochem.org.


LICENSE

Copyright (c) 2017 Raul Mera <rmera{at}usachDOTcl>


This program, including its documentation,
is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 2.0 as
published by the Free Software Foundation.

This program and its documentation is distributed in the hope that
it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General
Public License along with this program.  If not, see
<http://www.gnu.org/licenses/>.



*/

package main

import (
	"fmt"
	"math"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"

	"github.com/rmera/gochem/traj/stf"
	//	"sort"
)

type Closer interface {
	Close()
}

//Returns a function that superimposes each frame of a trajectory to the reference structure using the given atom,
//writing a DCD trajectory witht he result
func Super(mol *chem.Molecule, args []string, indexes []int) (func(coord *v3.Matrix) []float64, Closer) {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	var err error
	var neededargs int = 2
	if indexes == nil {
		if argslen < 1 {
			panic("Super: Got neither a selection nor a set of indexes!")
		}
		indexes, err = sel2atoms(mol, args[0])
		if err != nil {
			panic("Super: sel2atoms:" + err.Error())
		}
	}
	wname := "superimposed.stf"
	if argslen >= neededargs {
		wname = args[argslen-1]
	}
	wtraj, err := stf.NewWriter(wname, mol.Len(), nil) //I can't close this, sorry :D
	ret := func(coord *v3.Matrix) []float64 {
		super, err := chem.Super(coord, mol.Coords[0], indexes, indexes)
		if err != nil {
			panic("super: " + err.Error())
		}
		wtraj.WNext(super)
		return []float64{0.0} //meaningless
	}
	return ret, wtraj
}

func sSuper(ctest, ctempla, tmp *v3.Matrix) (float64, error) {
	if ctest.NVecs() != ctempla.NVecs() || tmp.NVecs() != ctest.NVecs() {
		return -1, fmt.Errorf("memRMSD: Ill formed matrices for memRMSD calculation")
	}
	tmp.Sub(ctest, ctempla)
	rmsd := tmp.Norm(2)
	return rmsd / math.Sqrt(float64(ctest.NVecs())), nil
}

//Returns a function that sums the position of each atom over over the trajectory and puts it into
//target, while saving the total frames read in the trajectory into N, so the average structure over
//the trajectory can be later obtained.
//I am a bit concerned about overflowing the float64s  in target
func Average(mol *chem.Molecule, target *v3.Matrix, N *int) func(coord *v3.Matrix) []float64 {
	ret := func(coord *v3.Matrix) []float64 {
		target.Add(coord, target)
		*N++
		return []float64{0.0} //meaningless
	}
	return ret
}
