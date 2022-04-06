package main

import (
	"log"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

//EasyShape takes a matrix of coordinates, a value for epsilon (a number close to zero, the closer, the more
//strict the orthogonality requriements are) and an (optative) masser and returns
//two shape indicators based on the elipsoid of inertia (or it massless equivalent)
//a linear and circular distortion indicators, and an error or nil (in that order).
//if you give a negative number as epsilon, the default (quite strict) will be used.
func FirstEig(coords *v3.Matrix, epsilon float64, mol ...chem.Masser) (*v3.Matrix, error) {
	var masses []float64
	var err error
	if len(mol) < 0 {
		masses = nil
	} else {
		masses, err = mol[0].Masses()
		if err != nil {
			masses = nil
			log.Println("FirstEig: Masses not found, will set them all to 1")
		}
	}
	moment, err := chem.MomentTensor(coords, masses)
	if err != nil {
		return nil, err
	}
	evecs, evals, err := v3.EigenWrap(moment, epsilon)
	if err != nil {
		return nil, err
	}
	log.Println(evals)

	return evecs.VecView(2), nil //The last and largest evec (I think!)

}

//Dimersangle is a specific function for SOD1 analysis. It will not go into the main
func DimersAngle(mol *chem.Molecule, args []string) func(coord *v3.Matrix) []float64 {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	epsilon := 0.0001
	if argslen < 4 {
		panic("DimersAngle: Not enough arguments, need 4!")
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
	//temporal vectors
	tp := v3.Zeros(1)
	tp2 := v3.Zeros(1)
	p := v3.Zeros(1)

	ret := func(coord *v3.Matrix) []float64 {
		vectors := make([]*v3.Matrix, 0, 2)
		planes := make([]*v3.Matrix, 0, 2)
		for i, v := range indexes {
			refmol.SomeAtoms(mol, v)
			refs[i].SomeVecs(coord, v)
			if i < 2 {
				tv, err := FirstEig(refs[i], epsilon, refmol)
				if err != nil {
					panic(err.Error)
				}
				vectors = append(vectors, tv)

			} else {
				tv, err := chem.BestPlane(refs[i], refmol)
				if err != nil {
					panic(err.Error)
				}
				planes = append(planes, tv)
			}
		}
		//We have the data, now comes the tricky part!
		ang := chem.Angle(planes[0], planes[1])
		ang *= chem.Rad2Deg
		if ang > 90 {
			planes[1].Scale(-1, planes[1]) //Should ensure that vectors are not "inverted"
		}
		tp.Add(planes[0], planes[1])
		p.Scale(0.5, tp)
		tp.Cross(p, vectors[0])
		tp2.Cross(p, tp)
		prv1 := chem.Projection(vectors[0], tp2)
		tp.Cross(p, vectors[1])
		tp2.Cross(p, tp)
		prv2 := chem.Projection(vectors[1], tp2)
		return []float64{chem.Angle(prv1, prv2) * chem.Rad2Deg}
	}
	return ret

}

//PlaneAngle function that will calculate the angle between the best planes for two selections, and do that for as many pairs of selections as requested.
func PlanesAngle(mol *chem.Molecule, args []string) func(coord *v3.Matrix) []float64 {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
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
