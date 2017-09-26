package main

import (
	"fmt"
	"github.com/rmera/gochem"
	"github.com/rmera/gochem/v3"
	"github.com/rmera/scu"
	"gonum.org/v1/gonum/mat"
	"sort"
	"strconv"
	"strings"
)

//ClosestN mol  "ref selection" "residue list" N

func ClosestN(mol *chem.Molecule, args []string) func(*v3.Matrix) []float64 {
	//	chem.FixGromacsPDB(mol)
	//	println("Hello from ClosestN!")
	argslen := len(args)
	if (argslen) < 3 {
		panic("ClosestN: Need at least two arguments")
	}
	N, err := strconv.Atoi(args[2])
	if err != nil {
		panic("ClosestN: " + err.Error())
	}
	refindexes, residues, com := resRankInput(mol, args)
	cutoff := 15.0 //We assume that N will be small enough that we won't go over 15 A. One can always change it, though.
	ret := func(coord *v3.Matrix) []float64 {
		retSlice := make([]float64, 0, N)
		ranked := distRank(coord, mol, refindexes, residues, cutoff, com)
		sort.Sort(MolDistList(ranked)) //possible bottleneck
		if len(ranked) < N {           //This shouldn't happen, but better have that.
			N = len(ranked)
		}
		for i := 0; i < N; i++ {
			retSlice = append(retSlice, ranked[i].Distance)
		}
		return retSlice
	}
	return ret
}

//returns a 1-member slice with the number of residues of the types specified of the selection found withint the given cutoff
func WithinCutoff(mol *chem.Molecule, args []string) func(*v3.Matrix) []float64 {
	//	chem.FixGromacsPDB(mol)
	//	println("Hello from ClosestN!")
	argslen := len(args)
	if (argslen) < 3 {
		panic("WithinCutoff: Need at least two arguments")
	}
	cutoff, err := strconv.ParseFloat(args[2], 64)
	if err != nil {
		panic("WithinCutoff: " + err.Error())
	}
	refindexes, residues, com := resRankInput(mol, args)
	ret := func(coord *v3.Matrix) []float64 {
		ranked := distRank(coord, mol, refindexes, residues, cutoff, com)
		//	fmt.Println(ranked)////////////////////////
		return []float64{float64(len(ranked))}
	}
	return ret
}

func resRankInput(mol *chem.Molecule, args []string) ([]int, []string, bool) {
	//	fmt.Println("Use: MDan distance sel1 sel2...")
	var err error
	//first we get the reference position from the first selection. If the selection has one atom
	//it is used, otherwise,
	refindex, err := sel2atoms(mol, args[0])
	if err != nil {
		panic("resRankInput: sel2atoms:" + err.Error())
	}
	residues := strings.Fields(args[1])
	com := false
	if scu.IsInString("--COM", args) { //the default will be to use each residues closest center of mass for calculating distances, but you can request to use the COM of the residue
		com = true
	}
	return refindex, residues, com
}

type resAndChain struct {
	ResID int
	Chain string
}

//returns all the residue numbers in mold covered by indexes
func allResIDandChains(mol chem.Atomer, indexes []int) []*resAndChain {
	ret := make([]*resAndChain, 0, 2)
	for _, i := range indexes {
		at := mol.Atom(i)
		if !repeated(at.MolID, at.Chain, ret) {
			ret = append(ret, &resAndChain{ResID: at.MolID, Chain: at.Chain})
		}
	}
	return ret
}

func repeated(id int, chain string, rac []*resAndChain) bool {
	for _, v := range rac {
		if v.Chain == chain && id == v.ResID {
			return true
		}
	}
	return false
}

func distRank(coord *v3.Matrix, mol *chem.Molecule, refindexes []int, residues []string, cutoff float64, com bool) []*MolDist {
	ranks := make([]*MolDist, 0, 30)
	resIDs := make([]*resAndChain, 0, 30)
	var molname string
	var id int
	var chain string
	var dist float64
	var err error
	water := v3.Zeros(3)
	var ref *v3.Matrix
	ownresIDs := allResIDandChains(mol, refindexes) //We could call this upstream and just get this numbers, but I suspect
	//it doesn't make much difference and the function already has lots of parameters. Something that can be improved if
	//the program is slow.
	if len(refindexes) == 1 { //no center of mass
		ref = coord.VecView(refindexes[0])
	} else {
		t1 := v3.Zeros(len(refindexes))
		ref, err = centerOfMass(coord, t1, mol, refindexes)
		if err != nil {
			panic("Distance: Func: " + err.Error())
		}
	}
	if cutoff < 0 {
		cutoff = 15 //if a negative cutoff is given we use 15 A as the default.
	}
	chunk := chem.NewTopology(0, 1)
	//	fmt.Println(chunk)////////////////
	for i := 0; i < mol.Len(); i++ {
		at := mol.Atom(i)
		molname = at.Molname
		id = at.MolID
		chain = at.Chain
		var test *v3.Matrix
		if scu.IsInString(molname, residues) && !repeated(id, chain, resIDs) && !repeated(id, chain, ownresIDs) { ////////////////////////////////
			//	fmt.Println(molname, residues)///////////////////////////////////////
			expectedreslen := 6
			indexes := make([]int, 1, expectedreslen)
			indexes[0] = i
			for j := i; j < i+80; j++ {
				at2 := mol.Atom(j)
				if at2.MolID != id {
					break
				}
				indexes = append(indexes, j)
				//	println("indexes res",id, j) ///////////////////////////////////
			}
			if len(indexes) == 3 {
				test = water
			} else {
				test = v3.Zeros(len(indexes))
			}
			test.SomeVecs(coord, indexes)
			if com {
				chunk.SomeAtoms(mol, indexes)
			}
			dist = molDist(test, ref, chunk, com)
			if dist <= cutoff {
				ranks = append(ranks, &MolDist{Distance: dist, MolID: id})
				resIDs = append(resIDs, &resAndChain{ResID: id, Chain: chain})
				//		break //////////////////////////////////////////////////////////////////////////
			}
		}
	}
	//	fmt.Println("NO WEI", len(resIDs), ranks) ///////////////////////////////////////////////////////////////////////
	return ranks

}

type MolDist struct {
	Distance float64
	MolID    int
}

func (M *MolDist) String() string {
	return fmt.Sprintf("D: %4.3f ID: %d", M.Distance, M.MolID)
}

type MolDistList []*MolDist

func (M MolDistList) Swap(i, j int) {
	M[i], M[j] = M[j], M[i]
}
func (M MolDistList) Less(i, j int) bool {
	return M[i].Distance < M[j].Distance
}
func (M MolDistList) Len() int {
	return len(M)
}

func molDist(test, ref *v3.Matrix, mol chem.Masser, com bool) float64 {
	temp := v3.Zeros(1)
	//	println("Call me!") /////////////
	if !com {
		var d1, dclosest float64
		var v1, vclosest *v3.Matrix
		vclosest = test.VecView(0)
		dclosest = dist(ref, vclosest, temp)
		for i := 1; i < test.NVecs(); i++ {
			v1 = test.VecView(i)
			d1 = dist(ref, v1, temp)
			//			fmt.Println(vclosest,v1,dclosest,d1)////////////////////////////
			if d1 < dclosest {
				vclosest = test.VecView(i)
				dclosest = d1
			}
		}
		return dclosest
	}
	mass, err := mol.Masses()
	if err != nil {
		panic("molDist: " + err.Error())
	}
	comass, err := centerOfMassII(test, mat.NewDense(test.NVecs(), 1, mass))
	if err != nil {
		panic("molDist: " + err.Error())
	}
	return dist(ref, comass, temp)
}

func dist(r, t, temp *v3.Matrix) float64 {
	temp.Sub(r, t)
	return temp.Norm(2)
}
