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
	"math"
	"os"

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

////use:  program [-skip=number -begin=number2] Task pdbname xtcname skip sel1 sel2 .... selN. Some tasks may require that N is odd that n is even.
func main() {
	//The skip options
	skip := flag.Int("skip", 0, "How many frames to skip between reads.")
	enforcemass := flag.Bool("enforcemass", false, "For tasks requiring atomic masses, exit the program if some masses are not available. Otherwise all masses are set to 1.0 if one or more values are not found.")

	begin := flag.Int("begin", 1, "The frame from where to start reading.")
	fixGromacs := flag.Bool("fixGMX", false, "Gromacs PDB numbering issue with more than 10000 residues will be fixed and a new PDB written")
	superTraj := flag.Bool("super", false, "No analysis is performed. Instead, the trajectory is superimposed to the reference structure")
	lovo := flag.Int("lovo", -1, "if >=0, uses LOVO to determine the residues used in a superposition. The number becomes the frames skipped during the LOVO calculation. See (and cite) 10.1371/journal.pone.0119264. This fag is only valid if the flag 'super' is set, and invalidates the 'tosuper' flag")
	lovolimit := flag.Float64("lovolimit", 1.0, "Only residues with a final RMSD less that this value (form a LOVO calculation), in A, will be considered for alignment. Only meaningful if the flag 'lovo' is set to >=0")
	tosuper := flag.String("tosuper", "", "The atoms to be used of the superposition, if that is to be performed")
	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "Usage:\n  %s: [flags] task geometry.pdb trajectory.xtc selection1 selection2 ... selectionN", os.Args[0])
		flag.PrintDefaults()
	}

	flag.Parse()
	args := flag.Args()
	//We do this first, as, if this task is given, we don't need any other parameter.
	if strings.Contains(strings.ToLower(args[0]), "selectionhelp") {
		fmt.Printf("The selections are defined in the following way: \"RESID1,RESID2,RESID3-RESID3+N,RESIDN CHAIN ATNAME1,ATNAME2\"\nRESID are residue numbers. They can be separated by commas or, to specify a range, with dashes: 12,13,120,125-128,145 That would select the residues 12,13,120,125,126,127,128,145\nCHAIN must be a chain identifier such as \"A\". If chain is \"ALL\", every chain will be used. \nATNAME is a PDB atom name such as CA (alpha carbon). Hydrogen names may vary with the forcefield employed. if ALL is given, as the first atom name, all atoms in the selected residues will be consiered.\n")
		os.Exit(0)
	}

	//	println("SKIP", *skip, *begin, args) ///////////////////////////
	var f func(*v3.Matrix) []float64
	mol, err := chem.PDBFileRead(args[1], false)
	if err != nil {
		panic(err.Error())
	}
	if *fixGromacs || strings.Contains(args[0], "ixGMX") {
		chem.FixGromacsPDB(mol)
		chem.PDBFileWrite("Fixed-"+args[1], mol.Coords[0], mol, nil)
		os.Exit(0) //no reason to keep going.
	}
	//If we don't find one or more masses, we just set them all to 1.0
	//unless you told us not to.
	_, err = mol.Masses()
	if err != nil && !*enforcemass {
		for i := 0; i < mol.Len(); i++ {
			at := mol.Atom(i)
			at.Mass = 1.0
		}
	}
	var superlist []int
	super := false
	if *tosuper != "" {
		superlist, err = sel2atoms(mol, *tosuper)
		if err != nil {
			panic("Wrong superposition list")
		}
		super = true
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
	if *superTraj {
		if *lovo >= 0 {
			fmt.Printf("LOVO alignment requested. Please cite: 10.1371/journal.pone.0119264\n")
			opt := align.DefaultOptions()
			opt.LessThanRMSD = *lovolimit
			name, chain := sel2nameandchain(args[3])
			opt.AtomNames = name
			opt.Chains = chain
			opt.Begin = *begin
			opt.Skip = *lovo
			fmt.Println(opt.AtomNames, opt.Chains, opt.Skip, opt.LessThanRMSD) /////////////////////
			fmt.Printf("Starting LOVO calculation. You might as well go for a coffee.\n")
			lovoret, err := align.LOVO(mol, mol.Coords[0], args[2], opt)
			fmt.Printf("LOVO calculation finished.\n")
			if err == nil {
				superlist = lovoret.Natoms //if it works, works
			} else {
				log.Printf("Couldn't obtain LOVO indexes for the superposition: %s", err.Error())
			}

		}
		var toclose Closer
		f, toclose = Super(mol, args[3:], superlist)
		defer toclose.Close()
		mdan(traj, mol.Coords[0], f, *skip, *begin, false, nil)
		return
	}

	//only used for the Average task
	var target *v3.Matrix
	var N int
	//used for rmsf
	var rmsf *rmsfstr

	task := strings.ToLower(args[0])
	if task == "distance" {
		f = Distance(mol, args[3:])
	} else if task == "rmsd" {
		f = RMSD(mol, args[3:])
	} else if task == "ramachandran" {
		f = Ramachandran(mol, args[3:])
	} else if task == "closestn" {
		f = ClosestN(mol, args[3:])
	} else if task == "rmsf" {
		f, rmsf = RMSF(mol, args[3:])
	} else if task == "withincutoff" {
		f = WithinCutoff(mol, args[3:])
	} else if task == "shape" {
		f = Shape(mol, args[3:])
	} else if task == "planesangle" {
		f = PlanesAngle(mol, args[3:])
	} else if task == "interbyres" {
		f = InterByRes(mol, args[3:])
	} else if task == "average" {
		target = v3.Zeros(mol.Len())
		f = Average(mol, target, &N)
	} else {
		fmt.Println("Args:", args)
		panic("Task parameter invalid or not present" + args[0])
	}
	mdan(traj, mol.Coords[0], f, *skip, *begin, super, superlist)
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
	//	println("Fuck you aaashooole")
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
		//		println("do I even run?") //////////////////////////////////
		limits := strings.Split(v, "-")
		//		println("limits",limits) ///////////////////////
		llimit, err := strconv.Atoi(limits[0])
		if err != nil {
			return nil, err
		}
		ulimit, err := strconv.Atoi(limits[1])
		if err != nil {
			return nil, err
		}
		//		println("llimit, ulimit", llimit, ulimit)
		for j := llimit; j <= ulimit; j++ {
			reslist = append(reslist, j)
		}
		//		println(reslist)           ////////////
		//		for _,v:=range(reslist){ //////////////////
		//			println(v)          //////////////////
		//		}               //////////////////////////////
	}
	//we get the chain:
	chain := fields[1]
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
		//	fmt.Println("selected", at.MolID, reslist[0],reslist[1],reslist[2]) /////////////////////////////////////////

		if (atnames[0] != "ALL" && !scu.IsInString(at.Name, atnames)) || (chain != "ALL" && at.Chain != chain) { //ALL makes all atoms in each residue to be included, or all chains.
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

/*****************RMSD family ***********/
//RMSD returns a function that will calculate the RMSD of as many selections as requested from a given set of coordinates against the coordinates
//in the mol object.
func RMSD(mol *chem.Molecule, args []string) func(coord *v3.Matrix) []float64 {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	if argslen < 1 {
		panic("RMSD: Not enough arguments, need at least one!")
	}
	indexes := make([][]int, 0, argslen)
	for _, v := range args {
		s, err := sel2atoms(mol, v)
		if err != nil {
			panic("RMSD: sel2atoms:" + err.Error())
		}
		indexes = append(indexes, s)
	}
	refs := make([]*v3.Matrix, 0, len(indexes))
	tests := make([]*v3.Matrix, 0, len(indexes))
	temps := make([]*v3.Matrix, 0, len(indexes))
	for _, v := range indexes {
		tr := v3.Zeros(len(v))
		tt := v3.Zeros(len(v))
		ttemp := v3.Zeros(len(v))
		tr.SomeVecs(mol.Coords[0], v) //the refs are already correctly filled
		refs = append(refs, tr)
		tests = append(tests, tt)
		temps = append(temps, ttemp)

	}
	ret := func(coord *v3.Matrix) []float64 {
		RMSDs := make([]float64, 0, len(indexes))
		for i, v := range indexes {
			tests[i].SomeVecs(coord, v) //the refs are already correctly filled
			rmsd, err := memRMSD(tests[i], refs[i], temps[i])
			if err != nil {
				panic("RMSD: " + err.Error())
			}
			RMSDs = append(RMSDs, rmsd)
		}
		return RMSDs
	}
	return ret
}

//smemRMSD calculates the RMSD between test and template, considering only the atoms
//present in the testlst and templalst for each object, respectively.
//It does not superimpose the objects.
//To save memory, it asks for the temporary matrix it needs to be supplied:
//tmp must be Nx3 where N is the number
//of elements in testlst and templalst
func memRMSD(ctest, ctempla, tmp *v3.Matrix) (float64, error) {
	if ctest.NVecs() != ctempla.NVecs() || tmp.NVecs() != ctest.NVecs() {
		return -1, fmt.Errorf("memRMSD: Ill formed matrices for memRMSD calculation")
	}
	tmp.Sub(ctest, ctempla)
	rmsd := tmp.Norm(2)
	return rmsd / math.Sqrt(float64(ctest.NVecs())), nil

}

/*******RMSF functions Family***************/

//RMSF returns a function that will calculate the RMSF of as many selections as requested from a given set of coordinates against the coordinates
//in the mol object.

type rmsfstr struct {
	data []string
}

func (r *rmsfstr) String() string {
	return strings.Join(r.data, "\n")
}

func RMSF(mol *chem.Molecule, args []string) (func(coord *v3.Matrix) []float64, *rmsfstr) {
	//	fmt.Println("Use: MDan RMSD sel1 sel2...")
	argslen := len(args)
	if argslen < 1 {
		panic("RMSF: Not enough arguments, need at least one!")
	}
	indexes := make([][]int, 0, argslen)
	for _, v := range args {
		s, err := sel2atoms(mol, v)
		if err != nil {
			panic("RMSF: sel2atoms:" + err.Error())
		}
		indexes = append(indexes, s)
	}
	frames := 0.0
	cm := make([]*v3.Matrix, 0, len(indexes))
	sqcm := make([]*v3.Matrix, 0, len(indexes))
	temp := make([]*v3.Matrix, 0, len(indexes))
	test := make([]*v3.Matrix, 0, len(indexes))
	//	stdevs:=make([][]float64,0,len(indexes))
	for _, v := range indexes {
		tr := v3.Zeros(len(v))
		tt := v3.Zeros(len(v))
		temptest := v3.Zeros(len(v))
		ttemp := v3.Zeros(len(v))
		tr.SomeVecs(mol.Coords[0], v) //the refs are already correctly filled
		cm = append(cm, tr)
		sqcm = append(sqcm, tt)
		temp = append(temp, ttemp)
		test = append(test, temptest)
		//stdevs=append(stdevs,make([]float64,0,len(v)))
	}
	rmsf := new(rmsfstr)
	ret := func(coord *v3.Matrix) []float64 {
		rmsf.data = make([]string, 0, len(indexes[0]))
		frames++
		numbers := 0
		//	output, _ := os.Create("RMSF.dat") //A bit crazy, but since I don't know when does the traj end, I have to write a "current" RMSF for each frame ( save for the first 2). I do it in the same file, which means that for each frame, the file gets overwritten.
		//		defer output.Close()
		for i, v := range indexes {
			test[i].SomeVecs(coord, v) //the refs are already correctly filled
			cm[i].Add(cm[i], test[i])
			temp[i].Dense.MulElem(test[i], test[i])
			sqcm[i].Add(temp[i], sqcm[i])
			if frames < 2 {
				continue
			}
			sqcm[i].Scale(1/frames, sqcm[i])
			cm[i].Scale(1/frames, cm[i])
			temp[i].MulElem(cm[i], cm[i])
			temp[i].Sub(sqcm[i], temp[i])

			//Since I don't know when do the frames stop, I need to every time get my accumulators back to the regular state.
			//Of course I could just multiply the new set of numbers to be added in each frame by 1/(frames-1), but I'll refrain from
			//getting cute until I know this works.
			sqcm[i].Scale(frames, sqcm[i])
			cm[i].Scale(frames, cm[i])
		}
		vecs := temp[0].NVecs()
		for j := 0; j < vecs; j++ {
			numbers++
			outstr := fmt.Sprintf("%7d ", numbers)
			for i := 0; i < len(temp); i++ {
				outstr = outstr + fmt.Sprintf("%8.3f ", math.Sqrt(temp[i].VecView(j).Norm(2)))
			}
			rmsf.data = append(rmsf.data, outstr)
		}
		return []float64{0.0} //Dummy output
	}
	return ret, rmsf
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
		if (lastread >= 0 && i < lastread+skip) || i < begin-1 { //not so nice check for this twice
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
