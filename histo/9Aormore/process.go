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
	upperlim := 9.0
	if len(os.Args) > 2 {
		upperlim = scu.MustParseFloat(os.Args[2], 64)
	}
	M1 := new(histo.Matrix)

	d1 := json.NewDecoder(f1)

	d1.Decode(M1)
	M1.NormalizeAll()
	//r, c := M1.Dims()
	divs := M1.CopyDividers()
	//dest := histo.NewMatrix(r, c, divs)
	//dest.Fill()
	//	histo.MatrixCombine(func(a, b, dest *histo.Data) { dest.Sub(a, b, true) }, M1, M2, dest)
	f := func(D *histo.Data) (float64, error) {
		d := D.View()
		if d == nil {
			return 0, nil
		}
		var ret float64
		for i, v := range d {
			if divs[i+1] >= upperlim {
				break
			}
			ret += v
		}
		return ret, nil
	}
	integrated, err := M1.FromAll(f)
	scu.QErr(err)
	fout, err := os.Create("histoproc.json")
	scu.QErr(err)
	defer fout.Close()
	j := json.NewEncoder(fout)
	j.Encode(integrated)
}
