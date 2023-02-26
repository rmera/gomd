package main

import (
	"encoding/json"
	"log"
	"os"
	"runtime"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/histo"
	v3 "github.com/rmera/gochem/v3"
	"github.com/rmera/scu"
)

func RamachandranHistogramMatrix(mol *chem.Molecule, args []string) (func(coord *v3.Matrix) []float64, *histo.Matrix) {
	var err error
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	step := 1.0 //degrees
	if argslen > 0 {
		step = scu.MustParseFloat(args[0])
	}
	first := 0
	last := 180
	ramaset, err := chem.RamaList(mol, "", []int{0, -1})
	scu.QErr(err)

	//this will be the matrix's rows and cols
	//now the dividers
	divs := make([]float64, 0, int(step)*last)
	for i := float64(first); i <= float64(last); i += step {
		divs = append(divs, i)
	}
	//the matrix will be Nx2, one col for each "rama residue" i.e
	matrix := histo.NewMatrix(len(ramaset), 2, divs)
	for i := 0; i < len(ramaset); i++ {
		for j := 0; j < 2; j++ {

			matrix.NewHisto(i, j, divs, nil, ramaset[i].MolID)
		}
	}
	ret := func(coord *v3.Matrix) []float64 {
		angles, err := chem.RamaCalc(coord, ramaset)
		if err != nil {
			panic("Ramachandran f: " + err.Error())
		}
		for i, v := range angles {
			matrix.AddData(i, 0, v[0])
			matrix.AddData(i, 1, v[1])
		}
		return []float64{0}
	}
	return ret, matrix
}

func prochisto(hist *histo.Matrix, name string) {
	fil, err := os.Create(name)
	js := json.NewEncoder(fil)
	js.Encode(hist)
	scu.QErr(err)
}

func COMDist(coord *v3.Matrix, indexes1, indexes2 []int, tmp1, tmp2 *v3.Matrix, tmpvecs [3]*v3.Matrix) float64 {
	var err error
	tmp1.SomeVecs(coord, indexes1)
	tmp2.SomeVecs(coord, indexes2)
	tmpvecs[0], err = chem.MassCenterMem(tmp1, tmp1, tmpvecs[0])
	scu.QErr(err)
	tmpvecs[1], err = chem.MassCenterMem(tmp2, tmp2, tmpvecs[1])
	scu.QErr(err)
	tmpvecs[2].Sub(tmpvecs[1], tmpvecs[0])
	//fmt.Println(tmp1, tmp2, tmpvecs, indexes1, indexes2, tmpvecs[2].Norm(2)) //////////////////
	return tmpvecs[2].Norm(2)
}

func COMDistPar(coord *v3.Matrix, indexes1, indexes2 []int) float64 {
	var err error
	tmp1 := v3.Zeros(len(indexes1))
	tmp2 := v3.Zeros(len(indexes2))
	tmpvecs := [3]*v3.Matrix{v3.Zeros(1), v3.Zeros(1), v3.Zeros(1)}
	tmp1.SomeVecs(coord, indexes1)
	tmp2.SomeVecs(coord, indexes2)
	tmpvecs[0], err = chem.MassCenterMem(tmp1, tmp1, tmpvecs[0])
	scu.QErr(err)
	tmpvecs[1], err = chem.MassCenterMem(tmp2, tmp2, tmpvecs[1])
	scu.QErr(err)
	tmpvecs[2].Sub(tmpvecs[1], tmpvecs[0])
	//fmt.Println(tmp1, tmp2, tmpvecs, indexes1, indexes2, tmpvecs[2].Norm(2)) //////////////////
	return tmpvecs[2].Norm(2)
}

/*****
Distance histogram functions

***/

//returns "BB" if at least an atom in mol has that name, or "CA" otherwise. Prints a warning to stderr.
func backboneName(mol chem.Atomer) string {
	log.Println("Name for the backbone atom not given, will use 'BB' (the name for the Martini backbone bead), if at least an atom with that name is present, or CA")
	var at *chem.Atom
	isbb := false
	name := "CA"
	for i := 0; i < mol.Len(); i++ {
		at = mol.Atom(i)
		if at.Name == "BB" {
			isbb = true
			break
		}
	}
	if isbb {
		name = "BB"
	}
	return name
}

func DistanceHistogramMatrix(mol *chem.Molecule, args []string) (func(coord *v3.Matrix) []float64, *histo.Matrix) {
	if len(args) < 1 {
		panic("GeneralHistogramMatrix needs at least 1 argument")
	}
	var f func(coord *v3.Matrix) []float64
	var m *histo.Matrix
	switch args[0] {
	case "centroid":
		f, m = CentroidDistanceHistogramMatrix(mol, args[1:])
	case "com":
		f, m = CentroidDistanceHistogramMatrix(mol, args[1:])
	default:
		f, m = BBDistanceHistogramMatrix(mol, args[1:])
	}
	return f, m
}

//CentroidDistanceHistogramMatrix returns a function that obtains a distance histogram between the centroids of  every pair of residues/submolecules in the mol molecule, along the trajectory. The histogram is also returned by DistanceHistogramMatrix but updated by the returned function as it goes over the trajectory.  Each residue needs to have one and only one atom with name args[0] (or "CA" by default)
func CentroidDistanceHistogramMatrix(mol *chem.Molecule, args []string) (func(coord *v3.Matrix) []float64, *histo.Matrix) {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	name := "CA"
	var at *chem.Atom
	if argslen < 1 {
		name = backboneName(mol)
	} else {
		name = args[0]
	}
	step := 0.1
	if argslen > 1 {
		step = scu.MustParseFloat(args[1])
	}
	first := 0
	last := 16
	if argslen > 2 {
		last = scu.MustAtoi(args[2]) //I'm lazy
	}
	res := make([]int, 0, mol.Atom(mol.Len()-1).MolID)
	chains := make([]string, 0, mol.Atom(mol.Len()-1).MolID)
	for i := 0; i < mol.Len(); i++ {
		at = mol.Atom(i)
		if at.Name != name {
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
	/*
		tmpmat := make([]*v3.Matrix, 0, len(indexes))
		for _, v := range indexes {
			tmpmat = append(tmpmat, v3.Zeros(len(v)))
		}
		tmpvecs := [3]*v3.Matrix{v3.Zeros(1), v3.Zeros(1), v3.Zeros(1)}
	*/
	//this will be the matrix's rows and cols
	Nres := len(indexes)
	//now the dividers
	divs := make([]float64, 0, int(step)*last)
	for i := float64(first); i <= float64(last); i += step {
		divs = append(divs, i)
	}
	matrix := histo.NewMatrix(Nres, Nres, divs)
	for i := 0; i < Nres; i++ {
		for j := 0; j < Nres; j++ {
			matrix.NewHisto(i, j, divs, nil)
		}
	}
	goruts := runtime.NumCPU() / 6
	wchans := make(chan bool, goruts)
	//working := make(chan bool, goruts)
	schans := make(chan bool, goruts)

	ret := func(coord *v3.Matrix) []float64 {

		for i, _ := range indexes {
			go func(in2 [][]int, curr int, ready, started chan bool) {
				//	println("pinche pendejo wey") ///////////////////
				started <- true
				tmpmat := make([]*v3.Matrix, 0, len(indexes))
				for _, v := range in2 {
					tmpmat = append(tmpmat, v3.Zeros(len(v)))
				}
				tmpvecs := [3]*v3.Matrix{v3.Zeros(1), v3.Zeros(1), v3.Zeros(1)}

				for j := 0; j < len(in2[:curr]); j++ {
					d := COMDist(coord, in2[curr], in2[j], tmpmat[curr], tmpmat[j], tmpvecs)
					matrix.AddDataCS(curr, j, d)
				}
				ready <- true
			}(indexes, i, wchans, schans)
			/*
				if len(wchans) >= 10 {
					for len(wchans) > 0 {
						<-wchans
					}
				}
			*/
			if len(schans) >= goruts {
				for len(schans) > 0 {
					<-schans
					<-wchans
				}
			}
		}
		defer func() {
			for len(schans) > 0 {
				<-schans
				<-wchans

			}
		}()

		return []float64{0}
	}

	return ret, matrix //the caller will have to deal with the final histogram.
	//It should be just one call to put it in a JSON file.
}

//BBDistanceHistogramMatrix returns a function that obtains a distance histogram from every pair of atoms in the mol molecule that have the desired name ("CA" by default), along the trajectory. The histogram is also returned by DistanceHistogramMatrix but updated by the returned function as it goes over the trajectory.
func BBDistanceHistogramMatrix(mol *chem.Molecule, args []string) (func(coord *v3.Matrix) []float64, *histo.Matrix) {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	name := "CA"
	if argslen < 1 {
		name = backboneName(mol)
	} else {
		name = args[0]
	}
	step := 0.1
	if argslen > 1 {
		step = scu.MustParseFloat(args[1])
	}
	first := 0
	last := 16
	if argslen > 2 {
		last = scu.MustAtoi(args[2]) //I'm lazy
	}

	//I need the number of residues
	backboneindex := make([]int, 0, 100)
	for i := 0; i < mol.Len(); i++ {
		at1 := mol.Atom(i)
		if at1.Name != name {
			continue
		}
		backboneindex = append(backboneindex, i)
	}
	//this will be the matrix's rows and cols
	Nres := len(backboneindex)
	//now the dividers
	divs := make([]float64, 0, int(step)*last)
	for i := float64(first); i <= float64(last); i += step {
		divs = append(divs, i)
	}
	matrix := histo.NewMatrix(Nres, Nres, divs)
	for i := 0; i < Nres; i++ {
		for j := 0; j < Nres; j++ {
			matrix.NewHisto(i, j, divs, nil)
		}
	}
	distvec := v3.Zeros(1)
	var v1, v2 *v3.Matrix
	ret := func(coord *v3.Matrix) []float64 {
		for i, v := range backboneindex {
			for j, w := range backboneindex[:i] { //inefficient
				v1 = coord.VecView(v)
				v2 = coord.VecView(w)
				distvec.Sub(v2, v1)
				//println(i, j) //////////////////
				matrix.AddData(i, j, distvec.Norm(2))
			}

		}
		return []float64{0}
	}
	return ret, matrix //the caller will have to deal with the final histogram.
	//It should be just one call to put it in a JSON file.
}
