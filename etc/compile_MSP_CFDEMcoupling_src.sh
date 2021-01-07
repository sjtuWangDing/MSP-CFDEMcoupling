#!/bin/bash

#===================================================================#
# Description: Compile MSP-CFDEMcoupling sources and will create all
#              lnInclude directories before compilation in order to
#              avoid missing headers in foreign libraries.
# Author: Wang Ding
#===================================================================#

source $MSP_CFDEM_TOOLS_DIR/compile_msp_cfdem_lib.sh

# create log direct
logDir="compile_log"
cd $MSP_CFDEM_ETC_DIR
mkdir -p $logDir

#================================================================================#
# must compile (but not clean) LIGGGHTS libraries
# compile src
#================================================================================#
compiledLibList="$MSP_CFDEM_ETC_DIR/package_msp_cfdem_lib_list.txt"
echo "Please provide the libraries to be compiled in the $compiledLibList file."

if [[ ! -f $compiledLibList ]]; then
  echo "$compiledLibList does not exist. Nothing will be done."
else
  # loop compiledLibList
  # 修改内部字段分隔符
  oldIFS=$IFS
  IFS=$'\n'
  for line in $(cat $compiledLibList); do
    if [[ "$line" = lib:* ]]; then
      # 从每一行的信息中获取到要编译的lib所在的目录名
      lineLen=$(expr length "$line")
      dirLen=$(expr $lineLen - 4)
      dirName=$(expr substr "$line" 5 $dirLen)
      path="$MSP_CFDEM_PROJECT_DIR/$dirName"
      if [ -d $path ]; then
        cd $path
        # linking include files to ./lnInclude
        wmakeLnInclude .
      fi
    fi
  done
  echo
  echo "Creation of lnInclude directories finished!"
  echo
  for line in $(cat $compiledLibList); do
    if [[ "$line" = lib:* ]]; then
      lineLen=$(expr length "$line")
      dirLen=$(expr $lineLen - 4)
      dirName=$(expr substr "$line" 5 $dirLen)
      path="$MSP_CFDEM_PROJECT_DIR/$dirName"
      if [ -d $path ]; then
        # compile lib
        logPath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"
        logFileName="log_compile_msp_cfdem_lib_$dirName"
        logFileContentHeader="$logFileName-$(date +%Y-%m-%d-%H:%M)"
        doClean=true
        compileMspCfdemLib $path $logPath $logFileName $logFileContentHeader $doClean
        if [ ${PIPESTATUS[0]} -ne 0 ]; then
          IFS=$oldIFS
          exit 1
        fi
      fi
    fi
  done
  IFS=$oldIFS
fi
