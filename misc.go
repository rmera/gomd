package main

import (
	"fmt"
	"log"
	"strings"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	"github.com/rmera/scu"
)

//CaptureVals is a function that captures the values of a file and writes a pdb file with
//the coordinates of the frame that has a value close to the reference value.
func CaptureVals(mol *chem.Molecule, args []string) func(coord *v3.Matrix) []float64 {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	if argslen < 2 {
		panic("CaptureVals: Not enough arguments, need at least one!")
	}
	valsfile, err := scu.NewMustReadFile(args[0])
	scu.QErr(err)
	refvalue := scu.MustParseFloat(args[1])
	tolerance := 0.5 //in A
	if argslen > 2 {
		tolerance = scu.MustParseFloat(args[2])
	}
	var data string
	ret := func(coord *v3.Matrix) []float64 {
		if data == "EOF" {
			return []float64{0.0}
		}
		data := valsfile.Next()
		if data == "EOF" {
			return []float64{0.0}
		}
		splitdata := strings.Fields(data)
		frame := splitdata[0]
		value := scu.MustParseFloat(splitdata[1])
		if !(value-tolerance <= refvalue) || !(value+tolerance >= refvalue) {
			return []float64{0.0}
		}
		name := fmt.Sprintf("%s_%5.3f", frame, value)
		name = strings.Replace(name, ".", "d", 1)
		name = name + ".pdb"
		chem.PDBFileWrite(name, coord, mol, nil)
		log.Println("Wrote", name)
		return []float64{0.0}
	}

	return ret
}
