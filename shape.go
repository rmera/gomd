package main

import (
	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

//PlaneAngle return function that will calculate the angle between the best plane for a selection in a test
// and a reference structure. The reference and test structures are assumed to be superimposed.
func PlaneAngle(mol *chem.Molecule, args []string) func(coord *v3.Matrix) []float64 {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	if argslen < 1 {
		panic("PlanesAngle: Need at least one argument!")
	}
	indexes := make([][]int, 0, argslen)
	for _, v := range args {
		s, err := sel2atoms(mol, v)
		if err != nil {
			panic("PlanesAngle: sel2atoms:" + err.Error())
		}
		indexes = append(indexes, s)
	}
	refmol1 := chem.NewTopology(0, 1)
	tests := make([]*v3.Matrix, 0, len(indexes)) //keeps temporary coordinates for plane calculations
	refs := make([]*v3.Matrix, 0, len(indexes))  //keeps the reference plane vectors
	for _, v := range indexes {
		tt := v3.Zeros(len(v))
		tests = append(tests, tt)
	}
	//we first obtain the vectors normal to the planes in the ref structure
	for i, v := range indexes {
		refmol1.SomeAtoms(mol, v)
		tests[i].SomeVecs(mol.Coords[0], v)
		p1, err := chem.BestPlane(tests[i], refmol1)
		if err != nil {
			panic("PlanesAngle: " + err.Error())
		}
		refs = append(refs, p1)

	}
	ret := func(coord *v3.Matrix) []float64 {
		pangles := make([]float64, 0, len(indexes))
		for i, v := range indexes {
			//for i, v := range indexes {
			refmol1.SomeAtoms(mol, v)
			tests[i].SomeVecs(coord, v)
			p2, err := chem.BestPlane(tests[i], refmol1)
			if err != nil {
				panic("PlanesAngle: " + err.Error())
			}
			ang := chem.Angle(refs[i], p2)
			ang *= chem.Rad2Deg
			if ang > 90 {
				ang = 180 - ang //We can't make this distinction, because the orientation of the vector normal to the moment-of-inertia plane can vary with small changes in the position of the atoms involved.
			}
			pangles = append(pangles, ang)

		}
		return pangles
	}
	return ret
}

//returns a function that will either obtain the angle between the planes for pairs of selections (if the "comp" option
//is given as the first option) or between the plane for a selection and the plane for the same selection in the reference
//structure (if the "comp" option is not given).
func PlanesAngle(mol *chem.Molecule, args []string) func(coord *v3.Matrix) []float64 {
	argslen := len(args)
	if argslen < 2 {
		panic("PlanesAngle: Need at least 1 arguments (3 if the 'comp' options is given)!")
	}
	if args[0] == "pairs" {
		return CompPlanesAngle(mol, args[1:])

	}
	return PlaneAngle(mol, args)

}

//CompPlanesAngle function that will calculate the angle between the best planes for two selections, and do that for as many pairs of selections as requested.
func CompPlanesAngle(mol *chem.Molecule, args []string) func(coord *v3.Matrix) []float64 {
	argslen := len(args)
	if argslen%2 != 0 || argslen < 2 {
		panic("PlanesAngle: Need a greater than 0 even number of arguments!")
	}
	indexes := make([][]int, 0, argslen)
	for _, v := range args {
		s, err := sel2atoms(mol, v)
		if err != nil {
			panic("PlanesAngle: sel2atoms:" + err.Error())
		}
		indexes = append(indexes, s)
	}
	refmol1 := chem.NewTopology(0, 1)
	refmol2 := chem.NewTopology(0, 1)
	refs := make([]*v3.Matrix, 0, len(indexes))

	for _, v := range indexes {
		tr := v3.Zeros(len(v))
		refs = append(refs, tr)
	}
	ret := func(coord *v3.Matrix) []float64 {
		pangles := make([]float64, 0, len(indexes))
		for i := 0; i < len(indexes)-1; i += 2 {
			//for i, v := range indexes {
			v := indexes[i]
			w := indexes[i+1]
			refmol1.SomeAtoms(mol, v)
			refmol2.SomeAtoms(mol, w)
			refs[i].SomeVecs(coord, v)   //the refs are already correctly filled
			refs[i+1].SomeVecs(coord, w) //the refs are already correctly filled
			p1, err := chem.BestPlane(refs[i], refmol1)
			if err != nil {
				panic("PlanesAngle: " + err.Error())
			}
			p2, err := chem.BestPlane(refs[i+1], refmol2)
			if err != nil {
				panic("PlanesAngle: " + err.Error())
			}
			ang := chem.Angle(p1, p2)
			ang *= chem.Rad2Deg
			if ang > 90 {
				ang = 180 - ang //We can't make this distinction, because the orientation of the vector normal to the moment-of-inertia plane can vary with small changes in the position of the atoms involved.
			}
			pangles = append(pangles, ang)

		}
		return pangles
	}
	return ret
}

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
