package main

import (
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	vor "github.com/rmera/gochem/voro"
	"github.com/rmera/scu"
)

//InterByRes returns a list with the residues involved in an interface, on each frame. As goMD functions needs to return floats,
//the list of residues will be in that format, but they will be integer numbers.
//this requires the "voro" package in goChem, which is not yet in master or devel.
//If both of the selections are from different chains, the residues from the second selection will be given as
//negative numbers, to distinguish them from residues from the other chain and the same residue ID..
func InterByRes(mol *chem.Molecule, args []string) func(coord *v3.Matrix) []float64 {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	if argslen < 2 {
		panic("InterByRes: I need exactly 2 selections and 1 optional float as arguments!")
	}
	indexA := wrapsel2atomsbyres(mol, args[0])
	indexB := wrapsel2atomsbyres(mol, args[1])
	indexT := make([]int, len(indexA), len(indexA)+len(indexB))
	copy(indexT, indexA)
	indexT = append(indexT, indexB...)
	//	fmt.Printf("lenA: %d lenB: %d lenT: %d\n", len(indexA), len(indexB), len(indexT)) ///////////
	//We will not be using the full molecule, but one where only atoms in indexT are considered, in order. Thus,
	//The first atom in indexA will now correspond to 0 in this new molecule, and so on. We adapt the indexes here.
	for i, _ := range indexA {
		indexA[i] = i
	}
	for i, _ := range indexB {
		indexB[i] = i + len(indexA)
	}
	//	fmt.Println(indexT) /////////////////
	nmol := chem.NewTopology(0, 1)
	nmol.SomeAtoms(mol, indexT)
	nmol.FillVdw()

	//	ncoord := v3.Zeros(len(indexT))
	cutoff := 4.0
	if argslen > 2 {
		var err error
		cutoff, err = strconv.ParseFloat(args[2], 64)
		if err != nil {
			panic("InterByRes: The third argument should be a number!: " + err.Error())
		}
	}
	chains1 := strings.Fields(args[0])[1]
	chains2 := strings.Fields(args[1])[1]
	c := strings.Contains

	//What follows is a stupid little trick to keep track of residues with the same ID but different chains
	//It's limited: It only supports one chain per selection (i.e. determine the interface between chains A and B, but not between a selection comprised of chains A and B and another comprised of chains A and C).
	//What I do is that, if there is one chain per selection and they are different, all the residue MolIDs from the second selection get multiplied by -1. It is easy afterwards to assign the negative MolIDs to the correct chain.
	chainmult := 1
	if !c(chains1, ",") && !c(chains2, ",") && chains1 != chains2 {
		chainmult = -1
	}

	ret := func(coord *v3.Matrix) []float64 {
		resinter := make([]int, 0, 30)

		ncoord := v3.Zeros(len(indexT))
		ncoord.SomeVecs(coord, indexT)
		planes := vor.GetPlanes(ncoord, nmol, cutoff, false)
		distv := v3.Zeros(1)
		for _, v := range indexA {
			at0 := nmol.Atom(v)
			//	if scu.IsInInt(at0.MolID, resinter) { //Residue already in the interface, no need to test this atom
			//		continue
			//	}
			for _, w := range indexB {
				at1 := nmol.Atom(w)
				v1 := ncoord.VecView(v)
				v2 := ncoord.VecView(w)
				distv.Sub(v2, v1)
				if distv.Norm(2) > 5 {
					continue
				}
				if scu.IsInInt(at0.MolID, resinter) && scu.IsInInt(at1.MolID*chainmult, resinter) { //Residues already in the interface
					continue
				}
				contact := planes.VdwContact(ncoord, nmol, v, w)
				if contact {
					resinter = append(resinter, at0.MolID)
					resinter = append(resinter, at1.MolID*chainmult)

				}

			}

		}

		return int2float(resinter)
	}
	return ret
}

//Plumbing

func wrapsel2atomsbyres(mol chem.Atomer, sel string) []int {
	in, err := sel2atoms(mol, sel)
	if err != nil {
		panic("InterByRes: sel2atoms: " + err.Error())
	}
	return in

}

//Returns a slice of floats with the same elements of the
//given slice of ints. Each repeated element is added only once!
func int2float(i []int) []float64 {
	f := make([]float64, 0, len(i))
	c := make([]int, 0, len(i))
	for _, v := range i {
		if scu.IsInInt(v, c) {
			continue
		}
		f = append(f, float64(v))
		c = append(c, v)
	}
	return f
}
