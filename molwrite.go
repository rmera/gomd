package main

import (
	"fmt"
	"log"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

func molprint(mol chem.Atomer, srange string) func(*v3.Matrix, []float64, int) {
	nilfunc := func(*v3.Matrix, []float64, int) {}
	prange := make([]float64, 2)
	if srange == "" {
		//Nothing wrong here, just that the user didn't request PDB printing.
		return nilfunc
	}
	splitrange := strings.Fields(srange)
	var err error
	if len(splitrange) != 2 {
		log.Printf("Couldn't parse range %s: %s. Will not print pdbs", srange, err.Error())
		return nilfunc
	}
	for i, v := range splitrange {
		prange[i], err = strconv.ParseFloat(v, 64)
		if err != nil {
			log.Printf("Couldn't parse range %s: %s. Will not print pdbs", srange, err.Error())
			return nilfunc
		}
	}
	f := func(coord *v3.Matrix, values []float64, index int) {
		if prange == nil {
			return
		}
		name := fmt.Sprintf("frame%d_", index)
		printed := false
		for i, v := range values {
			if v >= prange[0] && v <= prange[1] {
				name += fmt.Sprintf("valnumber%d_value%5.2f_", i+1, v)
				printed = true
			}
		}
		if printed {
			name += ".pdb"
			chem.PDBFileWrite(name, coord, mol, nil)

		}
	}
	return f
}
