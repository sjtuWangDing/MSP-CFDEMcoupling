#!/bin/bash

#===================================================================#
# Description: Compile LIGGGHTS libraries
# Author: Wang Ding
#===================================================================#

#- get time
timeNow="$(date +"%Y-%m-%d-%H:%M")"

#- create log direct
logDir="compile_log"
cd $MSP_CFDEM_ETC_DIR
mkdir -p $logDir

logPath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"

#================================================================================#
# compile src
#================================================================================#
whitelist="$MSP_CFDEM_ETC_DIR/package_undo_liggghts_list.txt"
echo ""
echo "deactivating all possible packages of LIGGGHTS now..."
echo "Please provide the packages to be compiled in the $whitelist file."
echo "Packages must be in: $LAMMPS_LIB_DIR."

if [[ ! -f $whitelist ]]; then
  echo "$whitelist does not exist in $MSP_CFDEM_ETC_DIR. Nothing will be done."
  nLines=0
  count=0
else
  nLines=`wc -l < $whitelist`
  count=0
fi

# resetting Makefile.package
cp $LIGGGHTS_SRC_DIR/Makefile.package.empty $LIGGGHTS_SRC_DIR/Makefile.package
cp $LIGGGHTS_SRC_DIR/Makefile.package.settings.empty $LIGGGHTS_SRC_DIR/Makefile.package.settings

while [ $count -lt $nLines ]
do
  let count++
  LINE=`head -n $count $whitelist | tail -1`
  if [[ "$LINE" == \#* ]]; then # comments
    continue
  elif [[ "$LINE" == "" ]]; then
    continue
  elif [[ "$LINE" == */dir ]]; then
    echo "will change path to:"
    LINE=$(echo "${LINE%????}") # ???? 代表去除 LINE 最后 4 个字符
    cd "$LIGGGHTS_SRC_DIR"
    echo $(pwd)
    logFileName="log_compile_$LINE""_lib"
    # delete log file
    if [[ -f $logFileName ]]; then
      rm $logFileName 2>&1 | tee -a $logPath/$logFileName
    fi
    make no-$LINE 2>&1 | tee -a $logPath/$logFileName
  fi
done

whitelist="$MSP_CFDEM_ETC_DIR/package_liggghts_list.txt"
echo ""
echo "activating packages of LIGGGHTS now..."
echo "Please provide the packages to be compiled in the $whitelist file."
echo "Packages must be in: $LAMMPS_LIB_DIR."

if [ ! -f $whitelist ]; then
  echo "$whitelist does not exist in $CWD. Nothing will be done."
  nLines=0
  count=0
else
  nLines=`wc -l < $whitelist`
  count=0
fi

while [ $count -lt $nLines ]
do
  let count++
  LINE=`head -n $count $whitelist | tail -1`
  if [[ "$LINE" == \#* ]]; then # comments
    continue
  elif [[ "$LINE" == "" ]]; then
    continue
  elif [[ "$LINE" == */dir ]]; then
    LINE=$(echo "${LINE%????}") # ???? 代表去除 LINE 最后 4 个字符
    logFileName="log_compile_$LINE""_lib"
    # assuming we need the poems lib if the package POEMS is activated
    if [[ "$LINE" == "POEMS" ]]; then
        echo "compile $LINE:"
        cd $LAMMPS_LIB_PATH/poems
        make -f Makefile.g++ clean 2>&1 | tee -a $logPath/$logFileName
        make -j $nProc -f Makefile.g++ lib 2>&1 | tee -a $logPath/$logFileName
    fi
  fi
done
