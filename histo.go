package main

import (
	"log"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/histo"
	v3 "github.com/rmera/gochem/v3"
	"github.com/rmera/scu"
)

//DistanceHistogramMatrix returns a function that obtains a distance histogram from every pair of atoms in the mol molecule that have the desired name ("CA" by default), along the trajectory. The histogram is also returned by DistanceHistogramMatrix but updated by the returned function as it goes over the trajectory.
func DistanceHistogramMatrix(mol *chem.Molecule, args []string) (func(coord *v3.Matrix) []float64, *histo.Matrix) {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	name := "CA"
	if argslen < 1 {
		log.Println("Name for the backbone atom not given, will use 'CA'")
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
		for i, v := range backboneindex[:len(backboneindex)-1] {
			for j, w := range backboneindex[i+1:] {
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
