/*
To the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche.

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
	//	"sort"
)

func Super(mol *chem.Molecule, args []string, target *[]*v3.Matrix) func(coord *v3.Matrix) []float64 {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	if argslen < 1 {
		panic("Super: Not enough arguments, need exactly one!")
	}
	indexes, err := sel2atoms(mol, args[0])
	if err != nil {
		panic("Super: sel2atoms:" + err.Error())
	}
	ret := func(coord *v3.Matrix) []float64 {
		super, err := chem.Super(coord, mol.Coords[0], indexes, indexes)
		if err != nil {
			panic("super: " + err.Error())
		}
		*target = append(*target, super)
		return []float64{0.0} //meaningless
	}
	return ret
}

//smemRMSD calculates the RMSD between test and template, considering only the atoms
//present in the testlst and templalst for each object, respectively.
//It does not superimpose the objects.
//To save memory, it asks for the temporary matrix it needs to be supplied:
//tmp must be Nx3 where N is the number
//of elements in testlst and templalst
func sSuper(ctest, ctempla, tmp *v3.Matrix) (float64, error) {
	if ctest.NVecs() != ctempla.NVecs() || tmp.NVecs() != ctest.NVecs() {
		return -1, fmt.Errorf("memRMSD: Ill formed matrices for memRMSD calculation")
	}
	tmp.Sub(ctest, ctempla)
	rmsd := tmp.Norm(2)
	return rmsd / math.Sqrt(float64(ctest.NVecs())), nil

}
