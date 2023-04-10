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
	"log"
	"math"
	"strings"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	"github.com/rmera/scu"
)

/**********RMSD family ***********/
//RMSD returns a function that will calculate the RMSD of as many selections as requested from a given set of coordinates against the coordinates
//in the mol object.
func RMSD(mol *chem.Molecule, args []string) func(coord *v3.Matrix) []float64 {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	if argslen < 1 {
		panic("RMSD: Not enough arguments, need at least one!")
	}
	indexes := make([][]int, 0, argslen)
	for _, v := range args {
		s, err := sel2atoms(mol, v)
		if err != nil {
			panic("RMSD: sel2atoms:" + err.Error())
		}
		indexes = append(indexes, s)
	}
	refs := make([]*v3.Matrix, 0, len(indexes))
	tests := make([]*v3.Matrix, 0, len(indexes))
	temps := make([]*v3.Matrix, 0, len(indexes))
	for _, v := range indexes {
		tr := v3.Zeros(len(v))
		tt := v3.Zeros(len(v))
		ttemp := v3.Zeros(len(v))
		tr.SomeVecs(mol.Coords[0], v) //the refs are already correctly filled
		refs = append(refs, tr)
		tests = append(tests, tt)
		temps = append(temps, ttemp)

	}
	ret := func(coord *v3.Matrix) []float64 {
		RMSDs := make([]float64, 0, len(indexes))
		for i, v := range indexes {
			tests[i].SomeVecs(coord, v) //the refs are already correctly filled
			rmsd, err := memRMSD(tests[i], refs[i], temps[i])
			if err != nil {
				panic("RMSD: " + err.Error())
			}
			RMSDs = append(RMSDs, rmsd)
		}
		return RMSDs
	}
	return ret
}

func residueIndexes(mol chem.Atomer, args []string) [][]int {
	argslen := len(args)
	name := "CA"
	allowedchains := "ALL"
	var at *chem.Atom
	if argslen > 1 {
		allowedchains = args[0]
	}
	if argslen < 2 {
		name = backboneName(mol)
	} else {
		name = args[1]
	}

	res := make([]int, 0, mol.Atom(mol.Len()-1).MolID)
	chains := make([]string, 0, mol.Atom(mol.Len()-1).MolID)
	for i := 0; i < mol.Len(); i++ {
		at = mol.Atom(i)
		if at.Name != name || (allowedchains != "ALL" && !strings.Contains(allowedchains, at.Chain)) {
			continue
		}
		res = append(res, at.MolID)
		chains = append(chains, at.Chain)
	}
	indexes := make([][]int, 0, len(res))
	for i, v := range res {
		in := chem.Molecules2Atoms(mol, []int{v}, []string{chains[i]})
		indexes = append(indexes, in)
	}
	return indexes
}

func COMCompare(test, ref *v3.Matrix, indexes1, indexes2 []int, tmpt, tmpr *v3.Matrix, tmpvecs [3]*v3.Matrix) float64 {
	var err error
	tmpt.SomeVecs(test, indexes1)
	tmpr.SomeVecs(ref, indexes2)
	tmpvecs[0], err = chem.MassCenterMem(tmpt, tmpt, tmpvecs[0])
	scu.QErr(err)
	tmpvecs[1], err = chem.MassCenterMem(tmpr, tmpr, tmpvecs[1])
	scu.QErr(err)
	tmpvecs[2].Sub(tmpvecs[1], tmpvecs[0])
	return tmpvecs[2].Norm(2)
}

//PerResidueRMSD returns a function that will return the RMSD for each residue
//it considers up to 3 arguments: The first one is a string with all the chains to be considered
// ex: "ABC". If not given, or the string "ALL" is given, all chains are considered. The second
//is the name of a "backbone" atom taht is found
//once and only once in every residue (if not given, the function will attempt to deduce it from the
//structure). The third, if given should be the string "com", which will cause the function to obtain,
//for each residue, the distance between the centroids of the residue in the reference and current structure.
func PerResidueRMSD(mol *chem.Molecule, args []string) func(coord *v3.Matrix) []float64 {
	indexes := residueIndexes(mol, args)
	var refs, tests []*v3.Matrix
	var com bool
	if len(args) > 2 && args[2] == "com" {
		com = true
		refs = make([]*v3.Matrix, 0, len(indexes))
		tests = make([]*v3.Matrix, 0, len(indexes))

	}
	tmps := make([]*v3.Matrix, 0, len(indexes))
	for _, v := range indexes {
		tmps = append(tmps, v3.Zeros(len(v)))
		if com {
			refs = append(refs, v3.Zeros(len(v)))
			tests = append(tests, v3.Zeros(len(v)))
		}
	}
	tmpvecs := [3]*v3.Matrix{v3.Zeros(1), v3.Zeros(1), v3.Zeros(1)}
	ret := func(coord *v3.Matrix) []float64 {
		var rmsd []float64
		for i, v := range indexes {
			if com {
				rmsd = append(rmsd, COMCompare(coord, mol.Coords[0], v, v, tests[i], refs[i], tmpvecs))
				continue
			}
			r, err := chem.MemRMSD(coord, mol.Coords[0], tmps[i], v, v)
			if err != nil {
				panic("PerResidueRMSD: " + err.Error())
			}

			rmsd = append(rmsd, r)
		}
		return rmsd
	}
	return ret
}

func PerAtomRMSD(mol *chem.Molecule, args []string) func(coord *v3.Matrix) []float64 {
	argslen := len(args)
	if argslen < 1 {
		panic("PerAtomRMSD: Not enough arguments, need one!")
	}
	if argslen > 1 {
		log.Printf("PerAtomRMSD: Too many arguments. Will ignore all but the first.")
	}
	indexes, err := sel2atoms(mol, args[0])
	if err != nil {
		panic("RMSD: sel2atoms:" + err.Error())
	}
	ref := v3.Zeros(len(indexes))
	test := v3.Zeros(len(indexes))
	temp := v3.Zeros(len(indexes))
	ref.SomeVecs(mol.Coords[0], indexes) //the refs are already correctly filled
	ret := func(coord *v3.Matrix) []float64 {
		var rmsd []float64
		var err error
		test.SomeVecs(coord, indexes) //the refs are already correctly filled
		rmsd, err = chem.MemPerAtomRMSD(test, ref, nil, nil, temp)
		if err != nil {
			panic("RMSD: " + err.Error())
		}
		return rmsd
	}
	return ret
}

//smemRMSD calculates the RMSD between test and template, considering only the atoms
//present in the testlst and templalst for each object, respectively.
//It does not superimpose the objects.
//To save memory, it asks for the temporary matrix it needs to be supplied:
//tmp must be Nx3 where N is the number
//of elements in testlst and templalst
func memRMSD(ctest, ctempla, tmp *v3.Matrix) (float64, error) {
	if ctest.NVecs() != ctempla.NVecs() || tmp.NVecs() != ctest.NVecs() {
		return -1, fmt.Errorf("memRMSD: Ill formed matrices for memRMSD calculation")
	}
	tmp.Sub(ctest, ctempla)
	rmsd := tmp.Norm(2)
	return rmsd / math.Sqrt(float64(ctest.NVecs())), nil

}

/*******RMSF functions Family***************/

//RMSF returns a function that will calculate the RMSF of as many selections as requested from a given set of coordinates against the coordinates
//in the mol object.

type rmsfstr struct {
	data []string
}

func (r *rmsfstr) String() string {
	return strings.Join(r.data, "\n")
}

func RMSF2(mol *chem.Molecule, args []string) (func(coord *v3.Matrix) []float64, *rmsfstr) {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	if argslen < 1 {
		panic("RMSF: Not enough arguments, need at least one!")
	}
	indexes := make([][]int, 0, argslen)
	for _, v := range args {
		s, err := sel2atoms(mol, v)
		if err != nil {
			panic("RMSF: sel2atoms:" + err.Error())
		}
		indexes = append(indexes, s)
	}
	frames := 0.0
	cm := make([]*v3.Matrix, 0, len(indexes))
	sqcm := make([]*v3.Matrix, 0, len(indexes))
	temp := make([]*v3.Matrix, 0, len(indexes))
	test := make([]*v3.Matrix, 0, len(indexes))
	//	stdevs:=make([][]float64,0,len(indexes))
	for _, v := range indexes {
		tr := v3.Zeros(len(v))
		tt := v3.Zeros(len(v))
		temptest := v3.Zeros(len(v))
		ttemp := v3.Zeros(len(v))
		tr.SomeVecs(mol.Coords[0], v) //the refs are already correctly filled
		cm = append(cm, tr)
		sqcm = append(sqcm, tt)
		temp = append(temp, ttemp)
		test = append(test, temptest)
		//stdevs=append(stdevs,make([]float64,0,len(v)))
	}
	rmsf := new(rmsfstr)
	ret := func(coord *v3.Matrix) []float64 {
		rmsf.data = make([]string, 0, len(indexes[0]))
		frames++
		numbers := 0
		//	output, _ := os.Create("RMSF.dat") //A bit crazy, but since I don't know when does the traj end, I have to write a "current" RMSF for each frame ( save for the first 2). I do it in the same file, which means that for each frame, the file gets overwritten.
		//		defer output.Close()
		for i, v := range indexes {
			test[i].SomeVecs(coord, v) //the refs are already correctly filled
			cm[i].Add(cm[i], test[i])
			temp[i].Dense.MulElem(test[i], test[i])
			sqcm[i].Add(temp[i], sqcm[i])
			if frames < 2 {
				continue
			}
			sqcm[i].Scale(1/frames, sqcm[i])
			cm[i].Scale(1/frames, cm[i])
			temp[i].MulElem(cm[i], cm[i])
			temp[i].Sub(sqcm[i], temp[i])

			//Since I don't know when do the frames stop, I need to every time get my accumulators back to the regular state.
			//Of course I could just multiply the new set of numbers to be added in each frame by 1/(frames-1), but I'll refrain from
			//getting cute until I know this works.
			sqcm[i].Scale(frames, sqcm[i])
			cm[i].Scale(frames, cm[i])
		}
		vecs := temp[0].NVecs()
		for j := 0; j < vecs; j++ {
			numbers++
			outstr := fmt.Sprintf("%7d ", numbers)
			for i := 0; i < len(temp); i++ {
				outstr = outstr + fmt.Sprintf("%8.3f ", math.Sqrt(temp[i].VecView(j).Norm(2)))
			}
			rmsf.data = append(rmsf.data, outstr)
		}
		return []float64{0.0} //Dummy output
	}
	return ret, rmsf
}

func RMSF(mol *chem.Molecule, args []string) (func(coord *v3.Matrix) []float64, *rmsfstr) {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	if argslen < 1 {
		panic("RMSF: Not enough arguments, need at least one!")
	}
	indexes := make([][]int, 0, argslen)
	for _, v := range args {
		s, err := sel2atoms(mol, v)
		if err != nil {
			panic("RMSF: sel2atoms:" + err.Error())
		}
		indexes = append(indexes, s)
	}
	frames := 0.0
	cm := make([]*v3.Matrix, 0, len(indexes))
	sqcm := make([]*v3.Matrix, 0, len(indexes))
	temp := make([]*v3.Matrix, 0, len(indexes))
	test := make([]*v3.Matrix, 0, len(indexes))
	//	stdevs:=make([][]float64,0,len(indexes))
	for _, v := range indexes {
		tr := v3.Zeros(len(v))
		tt := v3.Zeros(len(v))
		temptest := v3.Zeros(len(v))
		ttemp := v3.Zeros(len(v))
		tr.SomeVecs(mol.Coords[0], v) //the refs are already correctly filled
		cm = append(cm, tr)
		sqcm = append(sqcm, tt)
		temp = append(temp, ttemp)
		test = append(test, temptest)
		//stdevs=append(stdevs,make([]float64,0,len(v)))
	}
	rmsf := new(rmsfstr)
	ret := func(coord *v3.Matrix) []float64 {
		rmsf.data = make([]string, 0, len(indexes[0]))
		frames++
		numbers := 0
		//	output, _ := os.Create("RMSF.dat") //A bit crazy, but since I don't know when does the traj end, I have to write a "current" RMSF for each frame ( save for the first 2). I do it in the same file, which means that for each frame, the file gets overwritten.
		//		defer output.Close()
		for i, v := range indexes {
			test[i].SomeVecs(coord, v) //the refs are already correctly filled
			cm[i].Add(cm[i], test[i])
			temp[i].Dense.MulElem(test[i], test[i])
			sqcm[i].Add(temp[i], sqcm[i])
			if frames < 2 {
				continue
			}
			sqcm[i].Scale(1/frames, sqcm[i])
			cm[i].Scale(1/frames, cm[i])
			temp[i].MulElem(cm[i], cm[i])
			temp[i].Sub(sqcm[i], temp[i])

			//Since I don't know when do the frames stop, I need to every time get my accumulators back to the regular state.
			//Of course I could just multiply the new set of numbers to be added in each frame by 1/(frames-1), but I'll refrain from
			//getting cute until I know this works.
			sqcm[i].Scale(frames, sqcm[i])
			cm[i].Scale(frames, cm[i])
		}
		vecs := temp[0].NVecs()
		for j := 0; j < vecs; j++ {
			numbers++
			outstr := fmt.Sprintf("%7d ", numbers)
			for i := 0; i < len(temp); i++ {
				outstr = outstr + fmt.Sprintf("%8.3f ", math.Sqrt(temp[i].VecView(j).Norm(2)))
			}
			rmsf.data = append(rmsf.data, outstr)
		}
		return []float64{0.0} //Dummy output
	}
	return ret, rmsf
}
