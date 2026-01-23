#!/bin/bash

#SBATCH -J LFV
#SBATCH -p gpu,cpu,high_cpu -x gpu-0-2,compute-0-1
#SBATCH -N 1
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --open-mode=append
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb 
#SBATCH --comment python
#SBATCH --time 03:00:00
#SBATCH --hint=compute_bound
export X509_USER_PROXY=$HOME/proxy/x509up

echo "Using proxy at: $X509_USER_PROXY"
/opt/ohpc/pub/voms/voms-proxy-info --all

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
#source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc12-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh
echo "python skimonefile.py -Y $year -C $ch -I $infile -O ${outpath}/${outfile} 2>&1 | tee ${logdir}/${outfile%%root}log"
python skimonefile.py -Y $year -C $ch -I $infile -O ${outpath}/${outfile} 2>&1 | tee ${logdir}/${outfile%%root}log
