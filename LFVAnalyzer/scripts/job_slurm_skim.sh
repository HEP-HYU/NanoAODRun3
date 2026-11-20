#!/bin/bash

#SBATCH -J LFV
#SBATCH -p high_cpu,cpu,gpu -x gpu-0-1
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

export XRD_CONNECTIONRETRY=10
export XRD_REQUESTTIMEOUT=300
export XRD_STREAMTIMEOUT=300
export XRD_REDIRECTLIMIT=16

cd $workdir
#source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc12-opt/setup.sh  # Run2 UL (CentOS7)
#source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh      # Run3 early (EL9)
source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc15-opt/setup.sh       # Run3 current (EL9)

outroot="${outpath}/${outfile}"
logfile="${logdir}/${outfile%%root}log"

echo "python scripts/skimonefile.py -Y $year -C $ch -I $infile -O ${outroot} 2>&1 | tee ${logfile}"
python scripts/skimonefile.py -Y "$year" -C "$ch" -I "$infile" -O "${outroot}" 2>&1 | tee "${logfile}"
rc=${PIPESTATUS[0]}

if [ "${rc}" -ne 0 ]; then
    rm -f "${outroot}"
    alt_infile=""
    case "${infile}" in
        root://xrootd-cms.infn.it/*)
            alt_infile="root://cmsxrootd.fnal.gov/${infile#root://xrootd-cms.infn.it/}"
            ;;
        root://cmsxrootd.fnal.gov/*)
            alt_infile="root://xrootd-cms.infn.it/${infile#root://cmsxrootd.fnal.gov/}"
            ;;
    esac

    if [ -n "${alt_infile}" ] && { grep -qE "GetBasket|TNetXNGFile|Xrd|ReadBuffer|Cannot open|Failed to open" "${logfile}" || ! grep -q "\[ERROR\]" "${logfile}"; }; then
        echo "[INFO] Retrying skim (attempt 2/2) with fallback redirector: ${alt_infile}" > "${logfile}"
        python scripts/skimonefile.py -Y "$year" -C "$ch" -I "$alt_infile" -O "${outroot}" 2>&1 | tee -a "${logfile}"
        rc=${PIPESTATUS[0]}
        if [ "${rc}" -ne 0 ]; then
            rm -f "${outroot}"
        fi
    fi
fi

exit "${rc}"
