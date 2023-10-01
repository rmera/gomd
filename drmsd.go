/*

rmsd.go, part of gomd

This program makes extensive use of the goChem Computational Chemistry library.
If you use this program, we kindly ask you support it by to citing the library as:

R. Mera-Adasme, G. Savasci and J. Pesonen, "goChem, a library for computational chemistry", http://www.gochem.org.


LICENSE

Copyright (c) 2017-2021 Raul Mera <rmera{at}usachDOTcl>


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

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

// This is the exact same function as dRMSD. I have no idea why did I copied.
// I _suppose_ I was planning to change it, and never did. Just ignore it.
// dRMSD2 returns a function to obtain the "inter-monomer" or "dimer" RMSD, for, of course
// a molecule/protein which is at least dimeric.
// it takes 3 selections, the first one is the superposition selection for the first monomer
// the second is the superposition selection for the second monomer.
// The third, is the selection for determining the RMSD of one of the monomers
// (say, all alpha-carbons in that monomer).
// The dRMSD measure attempts to determine the inter-monomer deformation/motion
// excluding the internal deformation of each monomer.
func dRMSD2(mol *chem.Molecule, args []string) func(coord *v3.Matrix) []float64 {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	if argslen < 3 {
		panic("dRMSD: Not enough arguments, need at least 3!")
	}
	indexes := make([][]int, 0, argslen)
	for _, v := range args[:3] {
		s, err := sel2atoms(mol, v)
		if err != nil {
			panic("dRMSD: sel2atoms:" + err.Error())
		}
		indexes = append(indexes, s)
	}

	lovo1 := new(chem.Molecule)
	lovo1.Copy(mol)
	lovo2 := new(chem.Molecule)
	lovo2.Copy(mol)

	refs := make([]*v3.Matrix, 0, len(indexes))
	tests := make([]*v3.Matrix, 0, len(indexes))
	temp := v3.Zeros(len(indexes[2]))
	for _, v := range indexes {
		tr := v3.Zeros(len(v))
		tt := v3.Zeros(len(v))
		tr.SomeVecs(lovo1.Coords[0], v) //the refs are already correctly filled
		refs = append(refs, tr)
		tests = append(tests, tt)

	}
	ret := func(coord *v3.Matrix) []float64 {
		//	lovo1.Coords[0].Copy(mol.Coords[0])
		//	lovo2.Coords[0].Copy(mol.Coords[0])
		var err error
		//lovoA
		_, err = memSuper(lovo1.Coords[0], coord, refs[0], tests[0], indexes[0], indexes[0])
		if err != nil {
			panic("super: " + err.Error())
		}
		//lovoB
		_, err = memSuper(lovo2.Coords[0], coord, refs[1], tests[1], indexes[1], indexes[1])
		if err != nil {
			panic("super: " + err.Error())
		}
		tests[2].SomeVecs(lovo1.Coords[0], indexes[2])
		refs[2].SomeVecs(lovo2.Coords[0], indexes[2])

		rmsd, err := memRMSD(tests[2], refs[2], temp)
		return []float64{rmsd}
	}
	return ret
}

// dRMSD returns a function to obtain the "inter-monomer" or "dimer" RMSD, for, of course
// a molecule/protein which is at least dimeric.
// it takes 3 selections, the first one is the superposition selection for the first monomer
// the second is the superposition selection for the second monomer.
// The third, is the selection for determining the RMSD of one of the monomers
// (say, all alpha-carbons in that monomer).
// The dRMSD measure attempts to determine the inter-monomer deformation/motion
// excluding the internal deformation of each monomer.
func dRMSD(mol *chem.Molecule, args []string) func(coord *v3.Matrix) []float64 {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	if argslen < 3 {
		panic("dRMSD: Not enough arguments, need at least 3!")
	}
	//this is the same as the other function!
	indexes := make([][]int, 0, argslen)
	for _, v := range args[:3] {
		s, err := sel2atoms(mol, v)
		if err != nil {
			panic("dRMSD: sel2atoms:" + err.Error())
		}
		indexes = append(indexes, s)
	}

	lovo1 := new(chem.Molecule)
	lovo1.Copy(mol)
	lovo2 := new(chem.Molecule)
	lovo2.Copy(mol)

	refs := make([]*v3.Matrix, 0, len(indexes))
	tests := make([]*v3.Matrix, 0, len(indexes))
	temp := v3.Zeros(len(indexes[2]))
	for _, v := range indexes {
		tr := v3.Zeros(len(v))
		tt := v3.Zeros(len(v))
		tr.SomeVecs(lovo1.Coords[0], v) //the refs are already correctly filled
		refs = append(refs, tr)
		tests = append(tests, tt)

	}
	ret := func(coord *v3.Matrix) []float64 {
		//	lovo1.Coords[0].Copy(mol.Coords[0])
		//	lovo2.Coords[0].Copy(mol.Coords[0])
		var err error
		//lovoA
		_, err = memSuper(lovo1.Coords[0], coord, refs[0], tests[0], indexes[0], indexes[0])
		if err != nil {
			panic("super: " + err.Error())
		}
		//lovoB
		_, err = memSuper(lovo2.Coords[0], coord, refs[1], tests[1], indexes[1], indexes[1])
		if err != nil {
			panic("super: " + err.Error())
		}
		tests[2].SomeVecs(lovo1.Coords[0], indexes[2])
		refs[2].SomeVecs(lovo2.Coords[0], indexes[2])

		rmsd, err := memRMSD(tests[2], refs[2], temp)
		return []float64{rmsd}
	}
	return ret
}

// Super determines the best rotation and translations to superimpose the coords in test
// considering only the atoms present in the slices of int slices indexes.
// The first indexes slices will be assumed to contain test indexes and the second, template indexes.
// If you give only one, it will be assumed to correspond to test, if test has more atoms than
// elements on the indexes set, or templa, otherwise. If no indexes are given, all atoms on each system
// will be superimposed. The number of atoms superimposed on both systems must be equal.
// Super modifies the test matrix, but template and indexes are not touched.
func memSuper(test, templa, ctest, ctempla *v3.Matrix, indexes ...[]int) (*v3.Matrix, error) {
	if len(indexes) == 0 || indexes[0] == nil || len(indexes[0]) == 0 { //If you put the date in the SECOND slice, you are just messing with me.
		ctest = test
		ctempla = templa
	} else if len(indexes) == 1 {
		if test.NVecs() > len(indexes[0]) {
			//	ctest = v3.Zeros(len(indexes[0]))
			ctest.SomeVecs(test, indexes[0])
			ctempla = templa
		} else if templa.NVecs() > len(indexes[0]) {
			//	ctempla = v3.Zeros(len(indexes[0]))
			ctempla.SomeVecs(templa, indexes[0])
		} else {
			return nil, fmt.Errorf("chem.Super: Indexes don't match molecules")
		}
	} else {
		//	ctest = v3.Zeros(len(indexes[0]))
		ctest.SomeVecs(test, indexes[0])
		//	ctempla = v3.Zeros(len(indexes[1]))
		ctempla.SomeVecs(templa, indexes[1])
	}

	if ctest.NVecs() != ctempla.NVecs() {
		return nil, fmt.Errorf("chem.Super: Ill formed coordinates for Superposition")
	}

	_, rotation, trans1, trans2, err1 := chem.RotatorTranslatorToSuper(ctest, ctempla)
	if err1 != nil {
		return nil, err1
	}
	test.AddVec(test, trans1)
	//	fmt.Println("test1",test, rotation) /////////////77
	test.Mul(test, rotation)
	//	fmt.Println("test2",test) ///////////
	test.AddVec(test, trans2)
	//	fmt.Println("test3",test) ///////
	return test, nil
}
