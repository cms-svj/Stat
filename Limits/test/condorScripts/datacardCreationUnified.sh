#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh

xrdcp -s root://cmseos.fnal.gov//store/user/cfallon/CMSSW_10_2_13.tgz .
tar -xf CMSSW_10_2_13.tgz
rm CMSSW_10_2_13.tgz
export SCRAM_ARCH=slc6_amd64_gcc700
cd CMSSW_10_2_13/src/
scramv1 b ProjectRename
eval `scramv1 runtime -sh`

echo "CMSSW: "$CMSSW_BASE

ls -la Stat/Limits/python
cd Stat/Limits/test
#mkdir Fisher
#mkdir plots
#mkdir Residuals
echo "Arguments passed to this script are:"
echo "Name of output directory : ${1}"
echo "Mode: ${2}"
echo "Doing Systematics. ${7}"
if [ ${7} == "N" ]
then
  s="-s"
else
  s=""
fi

echo "Signal Parameters: ${3} ${4} ${5} ${6}"
cmd="python createDatacardsUnified.py -d ${1} -m ${2} -t -Z ${3} -D ${4} -R ${5} -A ${6} ${s}"


echo "combine commands:"
echo ${cmd}
echo ${cmd} >/dev/stderr

$cmd

SVJ_NAME="SVJ_mZprime${3}_mDark${4}_rinv${5}_alpha${6}"

# export items to EOS
echo "List all root files = "
ls *.root
echo "List all files"
ls 
echo "*******************************************"
EOSDIR=/store/user/cfallon/datacards_testSyst/${1}
OUTDIR=root://cmseos.fnal.gov/${EOSDIR}
echo "xrdcp output for condor"
for FILE in *.root *.pdf *.txt #Residuals/*.pdf plots/*.pdf Fisher/*.txt ${SVJ_NAME}/*.txt
do
  echo "xrdcp -f ${FILE} ${OUTDIR}${FILE}"
  xrdcp -f ${FILE} ${OUTDIR}${FILE} 2>&1
  rm ${FILE}
done
#rm -rf Residuals
#rm -rf plots
#rm -rf Fisher
#rm -rf ${SVJ_NAME}



cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_10_2_13
