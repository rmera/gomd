package main

import (
	"encoding/json"
	"fmt"
	"os"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/histo"
	"github.com/rmera/scu"
)

func main() {
	f1, err := os.Open(os.Args[1])
	scu.QErr(err)
	defer f1.Close()
	f2, err := os.Open(os.Args[2])
	scu.QErr(err)
	defer f2.Close()
	MolIDs1, err := ResList2MolID2(os.Args[3])
	scu.QErr(err)
	MolIDs2, err := ResList2MolID2(os.Args[4])
	scu.QErr(err)
	M1 := new(histo.Matrix)
	M2 := new(histo.Matrix)
	d1 := json.NewDecoder(f1)
	d2 := json.NewDecoder(f2)
	d1.Decode(M1)
	d2.Decode(M2)
	M1.NormalizeAll()
	M2.NormalizeAll()
	r, c := M2.Dims()
	divs := M2.CopyDividers()
	dest := histo.NewMatrix(r, c, divs)
	dest.Fill()
	histo.MatrixCombine(func(a, b, dest *histo.Data) { dest.Sub(a, b, true) }, M1, M2, dest)
	f := func(D *histo.Data) (float64, error) {
		d := D.View()
		if d == nil {
			return 0, nil
		}
		var ret float64
		for _, v := range d {
			ret += v
		}
		return ret, nil
	}
	integrated, err := dest.FromAll(f)
	scu.QErr(err)
	iii := scu.IsInInt
	ndists := 0
	accdist := 0.0
	for i, v := range integrated {
		for j, w := range v {
			if (iii(i, MolIDs1) && iii(j, MolIDs2)) || (iii(i, MolIDs2) && iii(j, MolIDs1)) {
				if w == 0 {
					continue
				}
				ndists++
				accdist += w
			}
		}
	}
	fmt.Println("Average absolute distance difference:", accdist/float64(ndists))
}

/**This is an attempt to clean this code to add it to goChem**/

func ResList2MolID2(residues string) ([]int, error) {
	reslistS := strings.Split(residues, ",")
	reslist := make([]int, 0, len(reslistS))
	for _, v := range reslistS {
		if !strings.Contains(v, "-") { //this is the simple case, this field is just a number that we parse and add to reslist
			val, err := strconv.Atoi(v)
			if err != nil {
				return nil, err
			}
			reslist = append(reslist, val)
			continue
		}
		//if there is a "-" in v, we have a bit more work to do.
		limits := strings.Split(v, "-")
		llimit, err := strconv.Atoi(limits[0])
		if err != nil {
			return nil, err
		}
		ulimit, err := strconv.Atoi(limits[1])
		if err != nil {
			return nil, err
		}
		for j := llimit; j <= ulimit; j++ {
			reslist = append(reslist, j)
		}
	}
	return reslist, nil
}

//Language for the selection, not very sophisticated:
// resIDs chain atoms
//resID is a list of residue IDs (integers) of one residue or a sequence of several, separated by commas and/or dashes. Numbers separated by dashed indicate that all numbers
//between the two (themselves included) are valid residue IDs that will be included.
//atoms is a list of atom names in the PDB/AMBER convention. Only atoms in this list will be included. ALL includes every atom.
//It the first character of atoms is an exclamation signe "!", all atoms but the ones in the list will be included.
//If a list of atom IDs (starting from 0), separated by spaces, and starting with the word "ATOMLIST" is given, the IDs will be returned. This is useful when processing a trajectory in xyz format which does not specify residues.
func Sel2Atoms(mol chem.Atomer, sel string) ([]int, error) {
	fields := strings.Fields(sel)
	//This will be used to parse a simple atomnombre list
	if fields[0] == "ATOMLIST" {
		indexes := make([]int, 0, len(fields[1:]))
		for _, i := range fields[1:] {
			index, err := strconv.Atoi(i)
			if err != nil {
				return nil, err
			}
			indexes = append(indexes, index)
		}
		return indexes, nil
	}
	//we start parsing the residue IDs
	reslist, err := ResList2MolID2(fields[0])
	if err != nil {
		return nil, err
	}

	//we get the chain:
	chain := strings.Split(fields[1], ",")
	atnames := strings.Split(fields[2], ",")
	ret := make([]int, 0, len(reslist)*4)
	for i := 0; i < mol.Len(); i++ {
		at := mol.Atom(i)
		if !scu.IsInInt(at.MolID, reslist) {
			continue
		}

		if (atnames[0] != "ALL" && !scu.IsInString(at.Name, atnames)) || (chain[0] != "ALL" && !scu.IsInString(at.Chain, chain)) { //ALL makes all atoms in each residue to be included, or all chains.
			continue
		}
		ret = append(ret, i)
		//		fmt.Println(ret)
	}
	return ret, nil
}
