#!/bin/bash

year=$1
ch=$2
filename=$3
infile="${filename}"
echo $infile
outpath=$4
outfile=$5
workdir=$6
logdir=$7

cd $workdir
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh
echo "python skimonefile.py -Y $year -C $ch -I $infile -O ${outpath}/${outfile} 2>&1 | tee ${logdir}/${outfile%%root}log"
python skimonefile.py -Y $year -C $ch -I $infile -O ${outpath}/${outfile} 2>&1 | tee ${logdir}/${outfile%%root}log
