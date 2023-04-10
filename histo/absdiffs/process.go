package main

import (
	"encoding/json"
	"os"

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
	fout, err := os.Create("histoproc.json")
	scu.QErr(err)
	defer fout.Close()
	j := json.NewEncoder(fout)
	j.Encode(integrated)
}
