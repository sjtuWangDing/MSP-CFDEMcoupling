#!/bin/bash:

#===================================================================#
# Description: Defining function to function to compile LIGGGHTS
# Author: Wang Ding
#===================================================================#

compileLIGGGHTS() {
  # define variables
  logPath="$1"
  logFileName="$2"
  logFileContentHeader="$3"
  doClean="$4"

  # clear old log file
  if [[ $logPath/$logFileName ]]; then
    rm $logPath/$logFileName
  fi

  # change path to $HOME/LIGGGHTS/LIGGGHTS-PUBLIC/src
  cd $LIGGGHTS_SRC_DIR

  # write content header to log file
  echo "\
/*---------------------------------------------------------------------------*\\
$logFileContentHeader
\*---------------------------------------------------------------------------*/
" 2>&1 | tee -a $logPath/$logFileName

  # write path
  echo "Log file path: $(pwd)" 2>&1 | tee -a $logPath/$logFileName
  # write empty line
  echo 2>&1 | tee -a $logPath/$logFileName

  # LIGGGHTS compilation flags
  if [[ $WM_NCOMPPROCS == "" ]]; then
    LIG_COMP_FLAGS=" "
  else
    LIG_COMP_FLAGS="-j $WM_NCOMPPROCS "
  fi

  # do clean and make
  if [[ $doClean == "false" ]]; then
    echo "Not cleaning LIGGGHTS"
  else
    echo "Cleaning LIGGGHTS: "
    make clean-$LIGGGHTS_MAKEFILE_NAME postfix=$LIGGGHTS_MAKEFILE_POSTFIX ${LIG_COMP_FLAGS} \
      2>&1 | tee -a $logPath/$logFileName
  fi

  echo "Compiling LIGGGHTS with ${LIG_COMP_FLAGS}:"
  make $LIGGGHTS_MAKEFILE_NAME postfix=$LIGGGHTS_MAKEFILE_POSTFIX ${LIG_COMP_FLAGS} \
    2>&1 | tee -a $logPath/$logFileName
  make makeshlib 2>&1 | tee -a $logPath/$logFileName
  make -f Makefile.shlib $LIGGGHTS_MAKEFILE_NAME postfix=$LIGGGHTS_MAKEFILE_POSTFIX ${LIG_COMP_FLAGS} \
    2>&1 | tee -a $logPath/$logFileName

  #   echo "compiling LIGGGHTS with ${LIG_COMP_FLAGS}"
  #   make $CFDEM_LIGGGHTS_MAKEFILE_NAME postfix=$CFDEM_LIGGGHTS_MAKEFILE_POSTFIX ${LIG_COMP_FLAGS} 2>&1 | tee -a $logpath/$logfileName
  #   make makeshlib 2>&1 | tee -a $logpath/$logfileName
  #   make -f Makefile.shlib $CFDEM_LIGGGHTS_MAKEFILE_NAME postfix=$CFDEM_LIGGGHTS_MAKEFILE_POSTFIX ${LIG_COMP_FLAGS} 2>&1 | tee -a $logpath/$logfileName

  if [ ${PIPESTATUS[0]} -ne 0 ]; then 
    return 1
  fi
}