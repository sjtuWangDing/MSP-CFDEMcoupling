#!/bin/bash

#===================================================================#
# Description: Execute CFD-DEM parallel run
# Author: Wang Ding
#===================================================================#

# source CFDEM env vars
. $HOME/.bashrc
. $MSP_CFDEM_bashrc

# set stack size to 100MB
defaultStackSize=`ulimit -s`
if [ $defaultStackSize -lt 102400 ]; then
  ulimit -s 102400
  echo "Set stack size to 100MB - done"
fi

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
if [ ! -d "$casePath/DEM/prev" ]; then
  mkdir "$casePath/DEM/prev"
fi

reInitDEM="false"
# if not exist restart file, execute DEM script to get restart file
if [[ ! -f "$casePath/DEM/prev/restart.init" || $reInitDEM = "true" ]]; then
  echo "execute DEM script to get restart file..."
  cd $casePath/DEM
  mpirun -np $numberOfProcsDEM liggghts < in.run
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

if [ $runOctave = "true" ]; then
  echo
  cd $casePath/CFD/octave
  octave --no-gui postproc.m
fi

if [ $postProc = "true" ]; then
  # get VTK data from liggghts dump file
  echo
  cd $casePath/DEM/post
  if [ ! -d "$casePath/DEM/post/coarse" ]; then
    mkdir "$casePath/DEM/post/coarse"
  fi
  if [ ! -d "$casePath/DEM/post/fine" ]; then
    mkdir "$casePath/DEM/post/fine"
  fi
  cp $casePath/DEM/post/dump_coarseParticles.run $casePath/DEM/post/coarse
  cp $casePath/DEM/post/dump_fineParticles.run $casePath/DEM/post/fine
  cd $casePath/DEM/post/coarse
  python2 $LPP_DIR/src/lpp.py dump_coarseParticles.run
  cd $casePath/DEM/post/fine
  python2 $LPP_DIR/src/lpp.py dump_fineParticles.run
fi

#- keep terminal open (if started in new terminal)
echo
echo "press Ctr+C kill process"
read
