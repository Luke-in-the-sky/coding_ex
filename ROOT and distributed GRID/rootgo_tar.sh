#!/bin/bash

### initialization
set echo
export PATH="${PATH}:/bin:/sbin:/usr/bin:/usr/local/bin:."

### dump some environment info
/bin/echo "PROCESS START"
/bin/echo ${HOSTNAME}
echo ${PATH}

### unpack source directory
/bin/tar -xzf rootgo.tgz
/bin/mv rootgo/* .
export WORKDIR=`pwd`

### copy data and root source from storage
/bin/echo "copy data file hits_0.raw"
lcg-cp --vo gilda lfn:/grid/gilda/lTriPadova/data2009/hits_0.raw file://`pwd`/hits_0.raw
if ! [ -e hits_0.raw ]; then exit 222; fi
/bin/echo "copy root source"
lcg-cp --vo gilda lfn:/grid/gilda/lTriPadova/data2009/rootsrc.tgz file://`pwd`/rootsrc.tgz
if ! [ -e rootsrc.tgz ]; then exit 223; fi
echo ""
echo " === Data retrieved "

### unpack root
/bin/ls -la
/bin/tar -xzf rootsrc.tgz

### install root
echo 
echo " === Installing root"
export ROOTSYS=${WORKDIR}/root
cd root
./configure linux
which g++
/usr/bin/gmake
echo
echo
echo
/bin/ls -la bin/root-config
/bin/ls -la bin/root
/bin/ls -la lib
if ! [ -e ./include/TFile.h ]; then exit 334; fi
if ! [ -e ${ROOTSYS}/include/TFile.h ]; then exit 335; fi

echo "fine root" >&2


### set ROOT environment
echo 
echo " === Setting environment "
cd ${WORKDIR}
echo ${ROOTSYS}
source envset
cd  ${ROOTSYS}/include
ls -lart
cd ${WORKDIR}
echo ${LIBLIST}
cp libRLDAP.so ${ROOTSYS}/lib

### compile and run
echo 
echo " === Compile: start" 
g++ -Wall -I${ROOTSYS}/include -o main main.cpp Track.cc WireDB.cc WireHit.cc Fitter.cc TrkTreeWrite.cc -L${ROOTSYS}/lib ${LIBLIST}
echo " === Executing"
./main hits_0.raw tracks-hits_0.root
