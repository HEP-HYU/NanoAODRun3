#!/bin/bash

#SBATCH -J LFV
#SBATCH -p gpu,cpu,high_cpu -x gpu-0-2
#SBATCH -N 1
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --open-mode=append
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb 
#SBATCH --comment python
#SBATCH --time 20:00:00
#SBATCH --hint=compute_bound

year=$1
ch=$2
indir=$3
outpath=$4
outfile=$5
workdir=$6
logdir=$7
syst=$8
mode=$9

cd $workdir
#source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc12-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh
echo "haha"

if [ "$#" -eq 8 ]; then
    if [[ "$indir" == *".root"* ]]; then
        echo "python processonefile.py -Y $year -C ${ch} -S ${syst} -I $indir -O ${outpath}/${outfile} 2>&1 | tee ${logdir}/${outfile%%root}log"
        python processonefile.py -Y $year -C ${ch} -S ${syst} -I $indir -O ${outpath}/${outfile} 2>&1 | tee ${logdir}/${outfile%%root}log
    else
        echo "python processonedataset.py -Y $year -C ${ch} -S ${syst} -I $indir -O ${outpath}/${outfile} 2>&1 | tee ${logdir}/${outfile%%root}log"
        python processonedataset.py -Y $year -C ${ch} -S ${syst} -I $indir -O ${outpath}/${outfile} 2>&1 | tee ${logdir}/${outfile%%root}log
    fi
fi
#FF for channel selction
if [ "$#" -eq 9 ]; then
    if [[ "$mode" == "--ff" ]]; then #FF application
        if [[ "$indir" == *".root"* ]]; then
            echo "python processonefile.py -Y $year -C ${ch} -S ${syst} -I $indir -O ${outpath}/${outfile} ${mode} 2>&1 | tee ${logdir}/${outfile%%root}log"
            python processonefile.py -Y $year -C ${ch} -S ${syst} -I $indir -O ${outpath}/${outfile} ${mode} 2>&1 | tee ${logdir}/${outfile%%root}log
        else
            echo "python processonedataset.py -Y $year -C ${ch} -S ${syst} -I $indir -O ${outpath}/${outfile} ${mode} 2>&1 | tee ${logdir}/${outfile%%root}log"
            python processonedataset.py -Y $year -C ${ch} -S ${syst} -I $indir -O ${outpath}/${outfile} ${mode} 2>&1 | tee ${logdir}/${outfile%%root}log
        fi
    else
        if [[ "$indir" == *".root"* ]]; then
            echo "python processonefile.py -Y $year -C ${ch} -S ${syst} -I $indir -O ${outpath}/${outfile} -M ${mode} 2>&1 | tee ${logdir}/${outfile%%root}log"
            python processonefile.py -Y $year -C ${ch} -S ${syst} -I $indir -O ${outpath}/${outfile} -M ${mode} 2>&1 | tee ${logdir}/${outfile%%root}log
        else
            echo "python processonedataset.py -Y $year -C ${ch} -S ${syst} -I $indir -O ${outpath}/${outfile} -M ${mode} 2>&1 | tee ${logdir}/${outfile%%root}log"
            python processonedataset.py -Y $year -C ${ch} -S ${syst} -I $indir -O ${outpath}/${outfile} -M ${mode} 2>&1 | tee ${logdir}/${outfile%%root}log
        fi
    fi
fi
