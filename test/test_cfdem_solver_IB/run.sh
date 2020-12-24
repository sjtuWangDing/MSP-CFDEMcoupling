#!/bin/bash

#===================================================================#
# Description: Execute CFD-DEM parallel run
# Author: Wang Ding
#===================================================================#

# source CFDEM env vars
. $HOME/.bashrc
. $MSP_CFDEM_bashrc

# include functions mspCfdemParallelRun
source $MSP_CFDEM_TOOLS_DIR/msp_cfdem_parallel_run.sh

# define variables for mpirun
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
solverDir="$MSP_CFDEM_BIN_DIR"
solverName="cfdemSolverIB-1.0"
numberOfProcs="4"
machineFileName="none" # yourMachinefileName
logPath=$casePath
logFileName="log_mpirun_$numberOfProcs_$solverName"

# block mesh
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
  echo "blockMesh: mesh was built before and using old mesh"
else
  echo "blockMesh: mesh needs to be built"
  cd $casePath/CFD
  blockMesh
fi

# add post dir to DEM
if [ ! -d "$casePath/DEM/post" ]; then
	mkdir "$casePath/DEM/post"
fi

# call function to run a parallel CFD-DEM case
mspCfdemParallelRun $casePath "$solverDir/$solverName" $numberOfProcs $machineFileName $logPath $logFileName

# define variables for post-process
runOctave="false"
postProc="false"
liggghtsDumpFileName="dump.liggghts_run"

if [ $runOctave = "true" ]; then
  echo
  cd $casePath/CFD/octave
  octave --no-gui postproc.m
fi

if [ $postProc = "true" ]; then
  # get VTK data from liggghts dump file
  echo
  cd $casePath/DEM/post
  python2 $LPP_DIR/src/lpp.py $liggghtsDumpFileName
fi

#- keep terminal open (if started in new terminal)
echo
echo "press Ctr+C kill process"
read
