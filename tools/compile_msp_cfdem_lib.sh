#!/bin/bash

# function to compile a msp-cfdem library
compileMspCfdemLib() {
  libPath="$1"
  logPath="$2"
  logFileName="$3"
  logFileContentHeader="$4"
  doClean="$5"

  # clear old log file
  if [ -e $logPath/$logFileName ]; then
    rm $logPath/$logFileName
  fi

  # write content header to log file
  echo "\
/*---------------------------------------------------------------------------*\\
$logFileContentHeader
\*---------------------------------------------------------------------------*/
" 2>&1 | tee -a $logPath/$logFileName
  # write path
  pwd 2>&1 | tee -a $logPath/$logFileName
  echo 2>&1 | tee -a $logPath/$logFileName

  cd "$libPath"
  if [ "true" = $doClean ]; then
    if [[ dev = $WM_PROJECT_VERSION || 3.0.* = $WM_PROJECT_VERSION || 4.* = $WM_PROJECT_VERSION || 5.* = $WM_PROJECT_VERSION ]]; then
      wrmdep 2>&1 | tee -a $logPath/$logFileName
    fi
    wclean 2>&1 | tee -a $logPath/$logFileName
  fi
  echo "Compiling a incompressible library."
  wmake 2>&1 | tee -a $logPath/$logFileName
  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    return 1
  fi
}
