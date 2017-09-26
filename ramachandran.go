package main

import (
	//	"fmt"
	"github.com/rmera/gochem"
	"github.com/rmera/gochem/v3"
	"math"
	"strconv"
	"strings"
)

//Ramachandran returns a Ramachandran-plotting function. For each residue in the selection
//it will plot 5 fields in each frame: Two for the Ramachandran plot and 3 for the corresponding RGB
//colors, so one can easily plot the Ramachandran evolution with simulation time.
//The arguments must be given as quoted strings containing a residue number and a chain. If RGB
//colors are to be included, the last argument must be "RGB" while the second-to-last argument
//must be the number of frames in the trajectory (needed to interpolate the colors).
func Ramachandran(mol *chem.Molecule, args []string) func(coord *v3.Matrix) []float64 {
	RGB := false
	var err error
	frames := -1 //the number of frames in the trajectory. -1 if unknown
	argslen := len(args)
	if argslen >= 3 && args[argslen-2] == "RGB" {
		frames, err = strconv.Atoi(args[argslen-1])
		if err != nil {
			panic("Ramachandran: " + err.Error())
		}
		args = args[0 : len(args)-2]
		argslen = argslen - 2
		RGB = true
	}
	ramasets := make([][]chem.RamaSet, 0, argslen)
	for _, v := range args {
		fields := strings.Fields(v)
		resid, err := strconv.Atoi(fields[0])
		if err != nil {
			panic("Ramachandran: " + err.Error())
		}
		set, err := chem.RamaList(mol, fields[1], []int{resid})
		if err != nil {
			panic("Ramachandran: " + err.Error())
		}
		ramasets = append(ramasets, set)
	}
	callnumber := 0
	var r, g, b float64
	ret := func(coord *v3.Matrix) []float64 {
		mult := 2
		if RGB {
			mult = 5
		}
		retSlice := make([]float64, 0, argslen*mult)
		for _, v := range ramasets {
			angles, err := chem.RamaCalc(coord, v)
			if err != nil {
				panic("Ramachandran f: " + err.Error())
			}
			retSlice = append(retSlice, angles[0][0])
			retSlice = append(retSlice, angles[0][1])
			if RGB {
				r, g, b = colors(callnumber, frames)

				retSlice = append(retSlice, r)
				retSlice = append(retSlice, g)
				retSlice = append(retSlice, b)
			}
			callnumber++
		}
		return retSlice
	}
	return ret
}

//Nice trick from icza on StackOverflow
//https://stackoverflow.com/questions/39544571/golang-round-to-nearest-0-05
//Note that id doesn't consider that when the rounded digit is five, numbers are promoted
//to the nearest even number. Still it's more than good enough for our purposes.
func round(x float64) float64 {
	return float64(int64(x + 0.5))
}

//These are the same goChem/chem functions. I will eventually make them
//exported in goChem and then replace these for calls to the original ones.
func colors(key, steps int) (r, g, b float64) {
	norm := 260.0 / float64(steps)
	hp := float64((float64(key) * norm) + 20.0)
	var h float64
	if hp < 55 {
		h = hp - 20.0
	} else {
		h = hp + 20.0
	}
	//	fmt.Println("HUE", h, hp)
	s := 1.0
	v := 1.0
	r, g, b = iHVS2RGB(h, v, s)
	return round(r), round(g), round(b) //slight change to the original to return rounded floats
}

//Here there is a change from the original gochem function to
//return float64 instead of uint8
//takes hue (0-360), v and s (0-1), returns r,g,b (0-255)
func iHVS2RGB(h, v, s float64) (float64, float64, float64) {
	var i, f, p, q, t float64
	var r, g, b float64
	maxcolor := 255.0
	conversion := maxcolor * v
	if s == 0.0 {
		return (conversion), (conversion), (conversion)
	}
	//conversion:=math.Sqrt(3*math.Pow(maxcolor,2))*v
	h = h / 60
	i = math.Floor(h)
	f = h - i
	p = v * (1 - s)
	q = v * (1 - s*f)
	t = v * (1 - s*(1-f))
	switch int(i) {
	case 0:
		r = v
		g = t
		b = p
	case 1:
		r = q
		g = v
		b = p
	case 2:
		r = p
		g = v
		b = t
	case 3:
		r = p
		g = q
		b = v
	case 4:
		r = t
		g = p
		b = v
	default: //case 5
		r = v
		g = p
		b = q
	}

	r = r * conversion
	g = g * conversion
	b = b * conversion
	return (r), (g), (b)
}
