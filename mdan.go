/*
To the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche.

goMD a little tool for the analysis of MD trajectories.


This program makes extensive use of the goChem Computational Chemistry library.
If you use this program, we kindly ask you support it by to citing the library as:

R. Mera-Adasme, G. Savasci and J. Pesonen, "goChem, a library for computational chemistry", http://www.gochem.org.


LICENSE

Copyright (c) 2017 Raul Mera <rmera{at}usachDOTcl>


This program, including its documentation,
is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 2.0 as
published by the Free Software Foundation.

This program and its documentation is distributed in the hope that
it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General
Public License along with this program.  If not, see
<http://www.gnu.org/licenses/>.



*/

package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"sort"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/align"
	"github.com/rmera/gochem/traj/amberold"
	"github.com/rmera/gochem/traj/dcd"
	"github.com/rmera/gochem/traj/stf"
	v3 "github.com/rmera/gochem/v3"
	"github.com/rmera/scu"
	"gonum.org/v1/gonum/mat"

	//	"sort"
	"strconv"
	"strings"
)

func formatLOVO(chain string, data []*align.MolIDandChain) string {
	first := true
	molids := make([][2]int, 1, 5)
	molids[0] = [2]int{0, 0}
	for _, v := range data {
		if chain != v.Chain() {
			continue
		}
		mid := v.MolID()
		if first {
			molids[0][0] = mid
			molids[0][1] = mid
			first = false
			continue
		}
		a := len(molids) - 1
		if mid == molids[a][1]+1 {
			molids[a][1] = mid
		} else {
			molids = append(molids, [2]int{mid, mid})
		}
	}
	ret := ""
	for _, v := range molids {
		if v[0] == v[1] {
			ret = fmt.Sprintf("%s,%d", ret, v[0])
		} else {
			ret = fmt.Sprintf("%s,%d-%d", ret, v[0], v[1])

		}

	}

	return fmt.Sprintf("%s %s CA", strings.TrimLeft(ret, ","), chain)

}

////use:  program [-skip=number -begin=number2] Task pdbname xtcname skip sel1 sel2 .... selN. Some tasks may require that N is odd that n is even.
func main() {
	pymol := flag.Bool("pymol", false, "Only for LOVO fit. Prints a command to create a pymol selection with the residues selected by the LOVO procedure.")
	nogmx := flag.Bool("gmx", false, "Only for LOVO fit. Do not print a command to create a Gromacs gmx make_ndx selection with the residues selected by the LOVO procedure. (this is printed by default)")

	skip := flag.Int("skip", 0, "How many frames to skip between reads.")
	enforcemass := flag.Bool("enforcemass", false, "For tasks requiring atomic masses, exit the program if some masses are not available. Otherwise all masses are set to 1.0 if one or more values are not found.")
	begin := flag.Int("begin", 1, "The frame from where to start reading.")
	lovo := flag.Int("lovo", -1, "if >=0, uses LOVO to determine the residues used in a superposition. The number becomes the frames skipped during the LOVO calculation. See (and cite) 10.1371/journal.pone.0119264. This flag invalidates the 'tosuper' flag")
	lovolimit := flag.Float64("lovolimit", 1.0, "Only residues with a final RMSD less that this value (form a LOVO calculation), in A, will be considered for alignment. Only meaningful if the flag 'lovo' is set to >=0")
	tosuper := flag.String("tosuper", "", "The atoms to be used of the superposition, if that is to be performed. Given as a regular goMD selection.")
	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "Usage:\n  %s: [flags] task geometry.pdb trajectory.xtc selection1 selection2 ... selectionN", os.Args[0])
		flag.PrintDefaults()
		fmt.Fprintf(flag.CommandLine.Output(), "\nAvailable tasks:  fixgmx, stop, distance, rmsd, peratomrmsd, ramachandran, closestn, rmsf, withincutoff, shape, planesangle, interbyres, average, super")
	}

	flag.Parse()
	args := flag.Args()
	//We do this first, as, if this task is given, we don't need any other parameter.
	if strings.Contains(strings.ToLower(args[0]), "selectionhelp") {
		fmt.Printf("The selections are defined in the following way: \"RESID1,RESID2,RESID3-RESID3+N,RESIDN CHAIN ATNAME1,ATNAME2\"\nRESID are residue numbers. They can be separated by commas or, to specify a range, with dashes: 12,13,120,125-128,145 That would select the residues 12,13,120,125,126,127,128,145\nCHAIN must be a chain identifier such as \"A\". If chain is \"ALL\", every chain will be used. \nATNAME is a PDB atom name such as CA (alpha carbon). Hydrogen names may vary with the forcefield employed. if ALL is given, as the first atom name, all atoms in the selected residues will be consiered.\n")
		os.Exit(0)
	}

	var f func(*v3.Matrix) []float64
	mol, err := chem.PDBFileRead(args[1], false)
	if err != nil {
		panic(err.Error())
	}
	//If we don't find one or more masses, we just set them all to 1.0
	//unless you told us not to. In the latter case, whatever task that
	//required massess wil crash.
	_, err = mol.Masses()
	if err != nil && !*enforcemass {
		for i := 0; i < mol.Len(); i++ {
			at := mol.Atom(i)
			at.Mass = 1.0
		}
	}

	var superlist []int
	if *tosuper != "" {
		superlist, err = sel2atoms(mol, *tosuper)
		if err != nil {
			panic("Wrong superposition list")
		}
	}
	var traj chem.Traj
	tmpf := strings.Split(args[2], ".")
	format := tmpf[len(tmpf)-1]
	switch strings.ToLower(format) {
	case "xtc":
		traj, err = OpenXTC(args[2]) //just a trick to offer
		if err != nil {
			panic(err.Error())
		}
	case "crd":
		traj, err = amberold.New(args[2], mol.Len(), false)
		if err != nil {
			panic(err.Error())
		}
	case "dcd":
		traj, err = dcd.New(args[2])
		if err != nil {
			panic(err.Error())
		}
	case "pdb":
		traj, err = chem.PDBFileRead(args[2], false)
		if err != nil {
			panic(err.Error())
		}
	case "xyz":
		traj, err = chem.XYZFileRead(args[2])
		if err != nil {
			panic(err.Error())
		}
	case "stz":
		traj, _, err = stf.New(args[2])
		if err != nil {
			panic(err.Error())
		}

	}
	var super bool = false

	//This is the only pre-processing we do.
	//You can request the LOVO procedure on its own (use the "stop" task, which does nothing) or use it to determine
	//the atoms for the "super" task.
	if *lovo >= 0 {
		fmt.Printf("#LOVO alignment requested. Please cite:\n# 10.1371/journal.pone.0119264\n# 10.1186/1471-2105-8-306\n")
		opt := align.DefaultOptions()
		opt.LessThanRMSD(*lovolimit)
		name, chain := sel2nameandchain(args[3])
		opt.AtomNames(name)
		opt.Chains(chain)
		opt.Begin(*begin)
		opt.Skip(*lovo)
		fmt.Println(opt.AtomNames(), opt.Chains(), opt.Skip(), opt.LessThanRMSD()) /////////////////////
		fmt.Printf("# Starting LOVO calculation. You might as well go for a coffee.\n")
		lovoret, err := align.LOVO(mol, mol.Coords[0], args[2], opt)
		fmt.Printf("# LOVO calculation finished.\n")
		if err == nil {
			superlist = lovoret.Natoms //if it works, works
		} else {
			log.Printf("Couldn't obtain LOVO indexes for the superposition: %s", err.Error())
		}
		sort.Sort(lovoret)
		//I'd like to have an option for VMD also, but don't know the selection 'language' there.
		fmt.Println("# LOVO atom indexes:", lovoret)
		if *pymol {
			fmt.Println("\n# PyMOL selection for LOVO-selected residues: ", lovoret.PyMOLSel())
		}
		if !*nogmx {
			fmt.Println("\n# gmx make_ndx selection LOVO-selected residues:\n ")
			for _, v := range chain {
				fmt.Println(lovoret.GMX(v))
			}
		}

		fmt.Println("# LOVO CA indexes in goMD selection format:")
		for _, v := range chain {
			fmt.Println(formatLOVO(v, lovoret.Nmols))
		}

	}

	//only used for the Average task
	var target *v3.Matrix
	var N int
	//used for rmsf
	var rmsf *rmsfstr
	task := strings.ToLower(args[0])
	switch task {
	case "fixgmx":
		chem.FixGromacsPDB(mol)
		chem.PDBFileWrite("Fixed-"+args[1], mol.Coords[0], mol, nil)
		return
	case "stop":
		return //nothing wrong, just something to use if you just want to do some of the pre-processing things
	case "distance":
		f = Distance(mol, args[3:])
	case "angle":
		f = Angle(mol, args[3:])
	case "rmsd":
		f = RMSD(mol, args[3:])
	case "peratomrmsd":
		f = PerAtomRMSD(mol, args[3:])
	case "ramachandran":
		f = Ramachandran(mol, args[3:])
	case "closestn":
		f = ClosestN(mol, args[3:])
	case "rmsf":
		f, rmsf = RMSF(mol, args[3:])
	case "withincutoff":
		f = WithinCutoff(mol, args[3:])
	case "shape":
		f = Shape(mol, args[3:])
	case "planesangle":
		f = PlanesAngle(mol, args[3:])
	case "interbyres":
		f = InterByRes(mol, args[3:])
	case "average":
		target = v3.Zeros(mol.Len())
		f = Average(mol, target, &N)
	case "super":
		var toclose Closer
		f, toclose = Super(mol, args[3:], superlist)
		defer toclose.Close()
		super = true
	default:
		fmt.Println("Args:", args)
		fmt.Println("Task parameter invalid or not present" + args[0])
		os.Exit(1)
	}
	mdan(traj, mol.Coords[0], f, *skip, *begin, super, superlist)

	//Extra post processing needed by some tasks

	if task == "average" && target != nil && N != 0 {
		cuoc := 1.0 / float64(N)
		target.Scale(cuoc, target)
		outname := strings.Replace(args[1], ".pdb", "_Average.pdb", -1)
		chem.PDBFileWrite(outname, target, mol, nil)

	}
	if task == "rmsf" {
		outp, err := os.Create("RMSF.dat")
		if err != nil {
			log.Printf("Couldn't open RMSF file for writing: %s", err.Error())
		}
		outp.WriteString(rmsf.String())
	}
}

/********General helper functions************/

func centerOfMass(coord, temp *v3.Matrix, mol *chem.Molecule, indexes []int) (*v3.Matrix, error) {
	top := chem.NewTopology(0, 1) //there is an goChem API change comming that will affect this call
	//	println("Fuck you aaashooole") //This is an Arnold reference. Sorry about the language :-)
	top.SomeAtoms(mol, indexes)
	temp.SomeVecs(coord, indexes)
	mass, err := top.Masses()
	if err != nil {
		return nil, err
	}
	//	println("Fuck YOU ashole!")
	//	fmt.Println(temp,indexes,mass) //////////////////
	ret, err := centerOfMassII(temp, mat.NewDense(len(indexes), 1, mass))
	//	fmt.Println("com!!",ret, temp, mass,indexes) ///////////////////////////
	return ret, err
}

//CenterOfMass returns the center of mass the atoms represented by the coordinates in geometry
//and the masses in mass, and an error. If mass is nil, it calculates the geometric center
func centerOfMassII(geometry *v3.Matrix, mass *mat.Dense) (*v3.Matrix, error) {
	if geometry == nil {
		return nil, fmt.Errorf("nil matrix to get the center of mass")
	}
	gr, _ := geometry.Dims()
	if mass == nil { //just obtain the geometric center
		tmp := ones(gr)
		mass = mat.NewDense(gr, 1, tmp) //gnOnes(gr, 1)
	}
	tmp2 := ones(gr)
	gnOnesvector := mat.NewDense(1, gr, tmp2)
	ref := v3.Zeros(gr)
	ref.ScaleByCol(geometry, mass)
	//	fmt.Println("ref", ref) ///////////////////
	ref2 := v3.Zeros(1)
	ref2.Mul(gnOnesvector, ref)
	ref2.Scale(1.0/mat.Sum(mass), ref2)
	return ref2, nil
}

//returns a flat64 slice of the size requested filed with ones
func ones(size int) []float64 {
	slice := make([]float64, size, size)
	for k, _ := range slice {
		slice[k] = 1.0
	}
	return slice
}

func sel2nameandchain(sel string) ([]string, []string) {
	var chains []string
	var names []string
	f := strings.Fields(sel)
	//first the chains
	//if chain is set to "ALL" we just return a nil slice, which means "consider all chains" for the LOVO function.
	if f[1] != "ALL" {
		chains = strings.Split(f[1], ",")
	}
	if f[2] == "ALL" {
		names = []string{"CA"} //yep, not having that "ALL" crap xD
	} else {
		names = strings.Split(f[2], ",")
	}
	return names, chains
}

//Language for the selection, not very sophisticated:
// resIDs chain atoms
//resID is a list of residue IDs (integers) of one residue or a sequence of several, separated by commas and/or dashes. Numbers separated by dashed indicate that all numbers
//between the two (themselves included) are valid residue IDs that will be included.
//atoms is a list of atom names in the PDB/AMBER convention. Only atoms in this list will be included. ALL includes every atom.
//It the first character of atoms is an exclamation signe "!", all atoms but the ones in the list will be included.
//If a list of atom IDs (starting from 0), separated by spaces, and starting with the word "ATOMLIST" is given, the IDs will be returned. This is useful when processing
//a trajectory in xyz format which does not specify residues and so on.
func sel2atoms(mol chem.Atomer, sel string) ([]int, error) {
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
	reslistS := strings.Split(fields[0], ",")
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
	//we get the chain:
	chain := strings.Split(fields[1], ",")
	//Now we go for the atomslist and chain
	/*
		negnames := false
		if strings.HasPrefix(fields[2], "!") {
			negnames = true
			fields[2] = fields[2][1:] //gotta check if I can do this and if the slicing is correct.
		}
	*/
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
	//We pick the atoms not present in ret. This is not working yet, so don't use it or included it in the documentation.
	/*
		if negnames {
			//ret are assumed to be sorted, which should be the case.
			ret2 := make([]int, 0, mol.Len()-len(ret))
			bs := 0 //where to begin the search in the ret string, i.e. the last index matched or tested.
			for i := 0; i < mol.Len(); i++ {
				if scu.IsInInt(i, ret[bs:]) {
					bs = i
					continue
				}
				bst := sort.SearchInts(ret[bs:], i)
				bs = bst + bs - 1 //the -1 is just in case I messed up to ensure we don't miss any index
				ret2 = append(ret2, i)

			}
			return ret2, nil
		}
	*/
	return ret, nil
}

/***************************
Now some of the Task functions (the rests are in other files)
****************************/

//********The Angle family functions**********//

func Angle(mol *chem.Molecule, args []string) func(*v3.Matrix) []float64 {
	//	fmt.Println("Use: MDan distance sel1 sel2...")
	argslen := len(args)
	if (argslen)%3 != 0 {
		panic("Distance: Always specity a number of selections that is divisible by 3")
	}
	com := make([]bool, 0, argslen/3) //should we use center of mass for this selection?
	indexes := make([][]int, 0, argslen)
	for i, v := range args {
		s, err := sel2atoms(mol, v)
		if err != nil {
			panic("Distance: sel2atoms:" + err.Error())
		}
		indexes = append(indexes, s)
		if i >= 2 && (i+1)%3 == 0 { //This ensure that we only check once for each angle "set" (i.e. every 3 selections)
			//	if len(indexes[i]) > 1 || len(indexes[i-1]) > 1 || len(indexes[i-2]) > 1 {
			if len(indexes[len(indexes)-1]) > 1 || len(indexes[len(indexes)-2]) > 1 || len(indexes[len(indexes)-3]) > 1 {
				com = append(com, true)
			} else {
				com = append(com, false)
			}
		}
	}
	//	fmt.Println(len(com), len(indexes), com, indexes) //////////////////////
	var vec1, vec2, vec3 *v3.Matrix

	//temporary angle vectors
	a1 := v3.Zeros(1)
	a2 := v3.Zeros(1)

	var err error
	//	fmt.Println(com) ////////////////////////////////////
	ret := func(coord *v3.Matrix) []float64 {
		angles := make([]float64, 0, len(indexes)/2)
		for i := 0; i < len(indexes); i = i + 3 { //we process them by triads
			if com[i/3] == false { //no center of mass indication, we assume taht each selection has one atom
				vec1 = coord.VecView(indexes[i][0])
				vec2 = coord.VecView(indexes[i+1][0])
				vec3 = coord.VecView(indexes[i+2][0])
			} else {
				//	println("get to the chooopaaaaa!")
				t1 := v3.Zeros(len(indexes[i]))
				t2 := v3.Zeros(len(indexes[i+1]))
				t3 := v3.Zeros(len(indexes[i+2]))

				vec1, err = centerOfMass(coord, t1, mol, indexes[i])
				if err != nil {
					panic("Distance: Func: " + err.Error())
				}
				vec2, err = centerOfMass(coord, t2, mol, indexes[i+1])
				if err != nil {
					panic("Distance: Func: " + err.Error())
				}
				vec3, err = centerOfMass(coord, t3, mol, indexes[i+1])
				if err != nil {
					panic("Distance: Func: " + err.Error())
				}

			}
			a1.Sub(vec1, vec2)
			a2.Sub(vec3, vec2)
			angles = append(angles, chem.Angle(a1, a2)*chem.Rad2Deg)
		}
		return angles
	}
	return ret
}

//********The Distance family functions**********//

func Distance(mol *chem.Molecule, args []string) func(*v3.Matrix) []float64 {
	//	fmt.Println("Use: MDan distance sel1 sel2...")
	argslen := len(args)
	if (argslen)%2 != 0 {
		panic("Distance: Always specity an even number of selections")
	}
	com := make([]bool, 0, argslen/2) //should we use center of mass for this selection?
	indexes := make([][]int, 0, argslen)
	for i, v := range args {
		s, err := sel2atoms(mol, v)
		if err != nil {
			panic("Distance: sel2atoms:" + err.Error())
		}
		indexes = append(indexes, s)
		if i >= 1 && (i+1)%2 == 0 { //we only need this value for
			if len(indexes[len(indexes)-1]) > 1 || len(indexes[len(indexes)-2]) > 1 { //tricky line, look here for errors.
				com = append(com, true)
			} else {
				com = append(com, false)
			}
		}
	}
	//	fmt.Println(len(com), len(indexes), com, indexes) //////////////////////
	var vec1, vec2 *v3.Matrix
	distvec := v3.Zeros(1) //the distance vector
	var err error
	//	fmt.Println(com) ////////////////////////////////////
	ret := func(coord *v3.Matrix) []float64 {
		distances := make([]float64, 0, len(indexes)/2)
		for i := 0; i < len(indexes); i = i + 2 { //we process them by pairs
			//		fmt.Println("com", i/2, com[i/2]) ///////////////////////////////////)
			if com[i/2] == false { //no center of mass indication, we assume taht each selection has one atom
				vec1 = coord.VecView(indexes[i][0])
				vec2 = coord.VecView(indexes[i+1][0])
			} else {
				//	println("get to the chooopaaaaa!")
				t1 := v3.Zeros(len(indexes[i]))
				t2 := v3.Zeros(len(indexes[i+1]))
				vec1, err = centerOfMass(coord, t1, mol, indexes[i])
				if err != nil {
					panic("Distance: Func: " + err.Error())
				}
				vec2, err = centerOfMass(coord, t2, mol, indexes[i+1])
				if err != nil {
					panic("Distance: Func: " + err.Error())
				}
				//			fmt.Println("Vectors",vec1, vec2)   ////////////////////
			}
			//			fmt.Println(vec2, vec1) /////////////////////////////////////////////
			distvec.Sub(vec2, vec1)
			distances = append(distances, distvec.Norm(2))
		}
		return distances
	}
	return ret
}

//Finally the "heart" of goMD, mdan.

//mdan takes a topology, a trajectory object and a function that must take a set of coordinates
//and a topology and returns a slice of floats. It applies the function to each snapshot of the trajectory.
//It then, for each snapshot, prints a line with the traj number as first field and the numbers in the returned
//slice as second to N fields, with the fields separated by spaces.
//the passed function should be a closure with everything necessary to obtain the desired data from each frame
//of the trajectory.
func mdan(traj chem.Traj, ref *v3.Matrix, f func(*v3.Matrix) []float64, skip, begin int, super bool, superlist []int) {
	var coords *v3.Matrix
	lastread := -1
	for i := 0; ; i++ { //infinite loop, we only break out of it by using "break"  //modified for profiling
		if lastread < 0 || (i >= lastread+skip && i >= begin-1) {
			coords = v3.Zeros(traj.Len())
		}
		err := traj.Next(coords) //Obtain the next frame of the trajectory.
		if err != nil {
			_, ok := err.(chem.LastFrameError)
			if ok || err.Error() == "EOF" {
				break //We processed all frames and are ready, not a real error.

			} else {
				panic(err.Error())
			}
		}
		if (lastread >= 0 && i < lastread+skip) || i < begin-1 { //not so nice having to check for this twice
			continue
		}
		if super {
			_, err := chem.Super(coords, ref, superlist, superlist)
			if err != nil {
				panic("Superposition failed! " + err.Error())
			}
		}
		lastread = i
		//The important part
		numbers := f(coords)
		fields := make([]string, len(numbers)+1)
		fields[0] = fmt.Sprintf("%7d", i+1)
		for j, v := range numbers {
			fields[j+1] = fmt.Sprintf("%8.3f", v)
		}
		str := strings.Join(fields, " ")
		fmt.Print(str + "\n")
		coords = nil // Not sure this works
	}

}
