#!/bin/bash

version=v0.7.1
XTC=1
xtcname=""
archi="Linux_x86-64"

rm -rf pack
mkdir pack

if [ "$XTC" -eq 1 ]
then
        echo "Compiling with xtc support"
        echo "The xtc-enabled binary will be named gomd"
        go build -tags "xtc"
        mv gomd pack/gomd
        xtcname="xtc"
else
        echo "Compiling without xtc support"
	xtcname=""
	go build
	mv gomd pack/gomd-nox
fi

go build



cp -R packagetools/xdrlib pack/
mkdir pack/plot
cp plot/ramachandran.py pack/plot/
cp plot/plots.py pack/plot/
cp plot/interfaces.py pack/plot/
cp plot/rdf2gomd.py pack/plot/
cp packagetools/INSTALL pack/
cp packagetools/gomd_config.sh pack/
cp gomd pack/
cp gomd pack/gomd-nox

cd pack
buildtime=$(date +%Y-%m-%d:%H:%M:%S)

tar -czvf ../gomd-$version-$archi-$xtcname-$buildtime.tgz *
cd ..



