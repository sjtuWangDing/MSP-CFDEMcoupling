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
solverName="mspCfdemSolverMix"
numberOfProcsDEM="16"
numberOfProcsCFD="8"
machineFileName="none" # yourMachinefileName
logPath=$casePath
logFileName="log_mpirun_$numberOfProcsCFD$solverName"

# add post dir to DEM
if [ ! -d "$casePath/DEM/post" ]; then
  mkdir "$casePath/DEM/post"
fi

# add post dir to DEM
if [ ! -d "$casePath/DEM/pre" ]; then
  mkdir "$casePath/DEM/pre"
fi

reInitDEM="false"
# if not exist restart file, execute DEM script to get restart file
if [[ ! -f "$casePath/DEM/pre/restart.liggghts_run" || $reInitDEM = "true" ]]; then
  echo "execute DEM script to get restart file..."
  cd $casePath/DEM
  mpirun -np $numberOfProcsDEM liggghts < in.liggghts_init
else
  echo "liggght: find restart file, no need to execute DEM script..."
fi

# block mesh
if [ -d "$casePath/CFD/constant/polyMesh/" ]; then
  echo "blockMesh: remove mesh..."
  rm -r $casePath/CFD/constant/polyMesh/
fi
echo "blockMesh: mesh will be built..."
cd $casePath/CFD
blockMesh

# call function to run a parallel CFD-DEM case
mspCfdemParallelRun $casePath "$solverDir/$solverName" $numberOfProcsCFD $machineFileName $logPath $logFileName

# define variables for post-process
runOctave="false"
postProc="true"
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
