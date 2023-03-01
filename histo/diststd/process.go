package main

import (
	"encoding/json"
	"fmt"
	"os"

	"github.com/rmera/gochem/histo"
	"github.com/rmera/scu"
)

func main() {
	f1, err := os.Open(os.Args[1])
	scu.QErr(err)
	defer f1.Close()
	M1 := new(histo.Matrix)
	d1 := json.NewDecoder(f1)
	d1.Decode(M1)
	M1.NormalizeAll()
	//r, c := M1.Dims()
	divs := M1.CopyDividers()
	fmt.Println("DIIIIIVS", divs) ////////////////////
	//dest := histo.NewMatrix(r, c, divs)
	//dest.Fill()
	//	histo.MatrixCombine(func(a, b, dest *histo.Data) { dest.Sub(a, b, true) }, M1, M2, dest)
	means, err := M1.FromAll(histo.Integrate(0, 6, divs))
	scu.QErr(err)
	fout, err := os.Create("means.json")
	scu.QErr(err)
	defer fout.Close()
	j := json.NewEncoder(fout)
	j.Encode(means)
	M1.UnNormalizeAll()
	stdev, err := M1.FromAll(histo.StdDev(divs))
	scu.QErr(err)
	M1.NormalizeAll()
	fout2, err := os.Create("stdevs.json")
	scu.QErr(err)
	defer fout2.Close()
	k := json.NewEncoder(fout2)
	k.Encode(stdev)
	median, err := M1.FromAll(histo.Quantile(0.5, divs))
	scu.QErr(err)
	fout3, err := os.Create("medians.json")
	scu.QErr(err)
	defer fout3.Close()
	l := json.NewEncoder(fout3)
	l.Encode(median)
	fmt.Println(stdev) ///////////////////////////////

}
