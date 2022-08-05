package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"runtime"
	"runtime/pprof"
	"strconv"
	"strings"

	"github.com/rmera/scu"
	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/stat"
)

//Remember kids, Stds are fun!
//Stds returns a slice with the standard deviation of each slice
//in vals
func Stds(vals [][]float64) []float64 {
	ret := make([]float64, 0, len(vals))
	for _, v := range vals {
		ret = append(ret, stat.StdDev(v, nil))
	}
	return ret

}

func mustfloat64(s string) float64 {
	f, err := strconv.ParseFloat(s, 64)
	if err != nil {
		panic(err.Error())
	}
	return f
}

func qerr(err error) {
	if err != nil {
		panic(err.Error())
	}
}

func makewindow(vals [][]float64, delay, blur int) [][]float64 {
	ret := make([][]float64, 0, len(vals))
	for _, v := range vals {
		ret = append(ret, make([]float64, 0, len(vals)))
		for j, _ := range v {
			if j >= delay {
				tmp := v[j-delay : j-delay+blur+1]
				ret[len(ret)-1] = append(ret[len(ret)-1], stat.Mean(tmp, nil))
			}
		}
	}
	return ret
}

var memprofile = flag.String("memprofile", "", "write memory profile to `file`")

func main() {
	cpus := flag.Int("cpus", 0, "CPU cores to be used. If 0 or less, all the available CPUs are used.")

	half := flag.Bool("half", false, "Obtain only half of the slopes/correlation matrix. Only to be used *without* delay. In that case, the matrix should be symmetrical, so half of the computer time is saved, and nothign is lost")

	nomax1 := flag.Bool("nomax1", false, "if slope >1 for the regression x vs y, the slope for regression y vs x is used. This option prevent that behavior. ")
	normalize := flag.Bool("n", false, "Normalize the RMSDs of each atom along the trajectory by dividing each value for the largest in the whole trajectory.")
	corr := flag.Bool("corr", false, "Obtain matrices of pearse-correlation coefficients")
	c := flag.Int("c", 306, "Columns in the gomd files, not counting the time (i.e first) column")
	//	filter := flag.Float64("filter", 0.0, "filter (set to zero) any value with absolute value smaller than the given")
	delay := flag.Int("delay", 0, "if >0, obtain a delayed-correlation index with a delay of the given number of frames")
	delayblur := flag.Int("delayblur", 0, "if --delay is used, for the 'delayed' data, use not the i-delay element, but the average between i-delay and i-delay+delayblur")
	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "Usage:\n  %s [flags] rsd.dat", os.Args[0])
		flag.PrintDefaults()
	}
	flag.Parse()
	args := flag.Args()
	if *cpus <= 0 {
		*cpus = runtime.NumCPU()
	}
	vals := make([][]float64, *c, *c)
	slopes := make([][]float64, *c, *c)
	for i := 0; i < *c; i++ {
		vals[i] = make([]float64, 0, 1000)
		slopes[i] = make([]float64, *c, *c)
	}
	fin, err := scu.NewMustReadFile(args[0])
	qerr(err)

	for i := fin.Next(); i != "EOF"; i = fin.Next() {
		l := strings.Fields(i)
		for j, _ := range vals {
			vals[j] = append(vals[j], mustfloat64(l[j+1])) //We skip the first "frame" column, so we count columns from 1
		}
		if *normalize {
			for i, v := range vals {
				floats.Scale(floats.Max(v), vals[i])
			}
		}
	}
	stds := Stds(vals)
	var disp [][]float64
	if *delay > 0 {
		disp = makewindow(vals, *delay, *delayblur)
	}
	var response func(x, y, w []float64, a bool) (float64, float64)
	correlation := func(x, y, w []float64, or bool) (float64, float64) {
		return 0, stat.Correlation(x, y, w)
	}
	response = stat.LinearRegression
	if *corr {
		response = correlation
	}
	chans := make([]chan bool, *cpus)
	for i := 0; i < *cpus; i++ {
		chans[i] = make(chan bool)
	}
	gorutines := 0
	for i, _ := range vals {
		go slicesslope(!*nomax1, *half, i, *delay, vals, disp, slopes, stds, response, chans[gorutines])
		gorutines++
		//The idea is that we pause the iteration when we have sent *cpu gorutines, and wait for them to return
		//The order doesn't quite matter, as they are all supposed to take give or take the same time.
		if gorutines == *cpus {
			for _, v := range chans {
				<-v
			}
			gorutines = 0
		}
	}
	//wait for any remaining gorutine to finish
	//Since the gorutines count is increased _after_ launching a gorutine, the index of the last channel passed to a gorutine will be gorutine-1
	for i := gorutines - 1; i >= 0; i-- {
		<-chans[i]
	}
	indicator := "slopes"
	if *corr {
		indicator = "corr"
	}
	foutname := fmt.Sprintf("%s_%s_del%0db%0d.json", strings.Replace(args[0], ".dat", "", -1), indicator, *delay, *delayblur)
	fout, err := os.Create(foutname)
	qerr(err)
	j := json.NewEncoder(fout)
	j.Encode(slopes)

	//mem profiling thing
	if *memprofile != "" {
		f, err := os.Create(*memprofile)
		if err != nil {
			log.Fatal("could not create memory profile: ", err)
		}
		defer f.Close() // error handling omitted for example
		runtime.GC()    // get up-to-date statistics
		if err := pprof.WriteHeapProfile(f); err != nil {
			log.Fatal("could not write memory profile: ", err)
		}
	}

}

//It is commonly safer to make anything that runs in a separate gorutine a proper, separate function (as opposed to a closure)
//That way you know that in only has access to the values you pass to it.
func slicesslope(max1, half bool, i, delay int, vals, disp, slopes [][]float64, stds []float64, response func(x, y, w []float64, a bool) (float64, float64), returnsignal chan bool) {
	var x, y []float64
	var beta float64
	for j, _ := range vals {
		if j <= i && half {
			continue
		}
		x = vals[i]
		y = vals[j]
		if delay > 0 {
			y = vals[j][delay:] //x is the delayed one, so I starts from a larger index to compensate
			x = disp[i]
		}
		// say, at the first position with delay 5 we will be comparing x[0] and y[4]
		stdv := stds[i] //these are wrong on the delayed case, but it doesn't matter. They are only here to save time.
		stdw := stds[j]
		//We always take the more "spread" set as x
		if stdv < stdw {
			x, y = y, x
		}
		_, beta = response(x, y, nil, false)
		//But we will revert if despite our choice, the slope is greater than 1 or smaller than -1
		if math.Abs(beta) > 1 && max1 {
			_, beta = response(y, x, nil, false)
		}
		slopes[i][j] = beta
	}
	//this is just a signal that says "I'm done". We don't care about the actual value.
	returnsignal <- true
}
