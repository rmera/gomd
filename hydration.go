package main

import (
	"fmt"
	"sort"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/solv"
	v3 "github.com/rmera/gochem/v3"
	"github.com/rmera/scu"
)

//RDF returns a function that will return
//the RDF or MDDF (depending on the number of reference atoms)
//for a set of atoms, up to the current frame.
func RDF(mol *chem.Molecule, args []string) func(*v3.Matrix) []float64 {
	//	chem.FixGromacsPDB(mol)
	//	println("Hello from ClosestN!")
	argslen := len(args)
	if (argslen) < 2 {
		panic("RDF: Need at least two arguments")
	}
	o := solv.DefaultOptions()
	o.Cpus(1)
	frames := 0
	refindexes, residues, com := resRankInput(mol, args)
	_ = com ////////gotta fix this
	totalsteps := int(o.End() / o.Step())
	acc := make([]float64, totalsteps)
	acctmp := make([]float64, totalsteps)
	distances := make([]float64, totalsteps)
	r := ""
	for i, _ := range distances {
		r = fmt.Sprintf("%s %3.5f", r, float64(i+1)*o.Step())
	}
	ret := func(coord *v3.Matrix) []float64 {
		mddf := solv.FrameUMolCRDF(coord, mol, refindexes, residues, o)
		frames++
		//	fmt.Println(mddf, frames) //////////
		for i, v := range mddf {
			acc[i] = acc[i] + v //this accumulates the values over frames  //(((frames - 1) * acc[i]) + v) / frames
			acctmp[i] = acc[i]  //this one gets resetted every frame
		}
		acctmp, _ = solv.MDFFromCDF(acctmp, frames, o.Step())
		fmt.Println(r) //The distances.
		return acctmp

	}
	return ret
}

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
	_ = com        ////////gotta fix this
	cutoff := 15.0 //We assume that N will be small enough that we won't go over 15 A. One can always change it, though.
	ret := func(coord *v3.Matrix) []float64 {
		retSlice := make([]float64, 0, N)
		ranked := distRank(coord, mol, refindexes, residues, cutoff)
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
	_ = com //should fix this
	ret := func(coord *v3.Matrix) []float64 {
		ranked := distRank(coord, mol, refindexes, residues, cutoff)
		return []float64{float64(len(ranked))}
	}
	return ret
}

func resRankInput(mol *chem.Molecule, args []string) ([]int, []string, bool) {
	var err error
	//first we get the reference position from the first selection. If the selection has one atom
	//it is used, otherwise,
	refindex, err := sel2atoms(mol, args[0])
	if err != nil {
		panic("resRankInput: sel2atoms:" + err.Error())
	}
	residues := strings.Fields(args[1])
	com := false
	if scu.IsInString("--COM", args) { //the default will be to use each residues closest atom for calculating distances, but you can request to use the COM of the residue
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

func distRank(coord *v3.Matrix, mol *chem.Molecule, refindexes []int, residues []string, cutoff float64) []*MolDist {
	ranks := make([]*MolDist, 0, 30)
	//	resIDs := make([]*resAndChain, 0, 30)
	var molname string
	var id, molid_skip int //The "skip" variables keep the residue just being read or discarded to avoid reading a residue twice.
	molid_skip = -1
	var chain, chain_skip string
	var distance float64
	water := v3.Zeros(3)
	var ref *v3.Matrix
	ownresIDs := allResIDandChains(mol, refindexes) //We could call this upstream and just get this numbers, but I suspect
	ref = v3.Zeros(len(refindexes))
	ref.SomeVecs(coord, refindexes)
	if cutoff < 0 {
		cutoff = 10 //if a negative cutoff is given we use 10 A as the default.
	}
	cutoffplus := cutoff + 2
	chunk := chem.NewTopology(0, 1)
	tmp := v3.Zeros(1)
	for i := 0; i < mol.Len(); i++ {
		at := mol.Atom(i)
		molname = at.Molname
		id = at.MolID
		chain = at.Chain
		var test *v3.Matrix
		//a little pre-screening for waters
		if scu.IsInString(molname, []string{"SOL", "WAT", "HOH"}) && dist(ref.VecView(0), coord.VecView(i), tmp) > cutoffplus {
			continue

		}
		if scu.IsInString(molname, residues) && !repeated(id, chain, ownresIDs) && (id != molid_skip || chain != chain_skip) {
			expectedreslen := 6
			indexes := make([]int, 1, expectedreslen)
			indexes[0] = i
			for j := i + 1; j < i+80; j++ {
				at2 := mol.Atom(j)
				if at2.MolID != id {
					break
				}
				indexes = append(indexes, j)
			}
			if len(indexes) == 3 {
				test = water
			} else {
				test = v3.Zeros(len(indexes))
			}
			test.SomeVecs(coord, indexes)
			distance = molDist(test, ref, chunk)
			if distance <= cutoff {
				ranks = append(ranks, &MolDist{Distance: distance, MolID: id})
				molid_skip = id
				chain_skip = chain
			}
		}
	}
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

func molDist(test, ref *v3.Matrix, mol chem.Masser) float64 {
	temp := v3.Zeros(1)
	var d1, dclosest float64
	var vt1, vr1 *v3.Matrix // vtclosest,vr1, vrclosest *v3.Matrix
	dclosest = dist(ref.VecView(0), test.VecView(0), temp)
	for i := 1; i < test.NVecs(); i++ { //This will fail if the test residues have only one atom.
		vt1 = test.VecView(i)
		for j := 0; j < ref.NVecs(); j++ {
			vr1 = ref.VecView(j)
			d1 = dist(vr1, vt1, temp)
			if d1 < dclosest {
				dclosest = d1
			}
		}
	}
	return dclosest
}

func dist(r, t, temp *v3.Matrix) float64 {
	temp.Sub(r, t)
	return temp.Norm(2)
}
