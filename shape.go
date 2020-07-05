package main

import (
	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

//Shape returns a function that will calculate shape indicators (planarity and elongation, returned in that order) on as many selections as requested
func Shape(mol *chem.Molecule, args []string) func(coord *v3.Matrix) []float64 {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	epsilon := 0.0001
	if argslen < 1 {
		panic("Shape: Not enough arguments, need at least one!")
	}
	indexes := make([][]int, 0, argslen)
	for _, v := range args {
		s, err := sel2atoms(mol, v)
		if err != nil {
			panic("Shape: sel2atoms:" + err.Error())
		}
		indexes = append(indexes, s)
	}
	refmol := chem.NewTopology(0, 1)
	refs := make([]*v3.Matrix, 0, len(indexes))

	for _, v := range indexes {
		tr := v3.Zeros(len(v))
		refs = append(refs, tr)
	}
	ret := func(coord *v3.Matrix) []float64 {
		shapes := make([]float64, 0, len(indexes))
		for i, v := range indexes {
			refmol.SomeAtoms(mol, v)
			refs[i].SomeVecs(coord, v) //the refs are already correctly filled
			lin, circ, err := chem.EasyShape(refs[i], epsilon, refmol)
			if err != nil {
				panic("Shape: " + err.Error())
			}
			shapes = append(shapes, circ)
			shapes = append(shapes, lin)
		}
		return shapes
	}
	return ret
}
