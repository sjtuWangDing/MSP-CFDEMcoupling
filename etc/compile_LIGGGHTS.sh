#!/bin/bash

#===================================================================#
# Description: Compile LIGGGHTS
# Author: Wang Ding
#===================================================================

#- include functions
source $MSP_CFDEM_TOOLS_DIR/compile_LIGGGHTS.sh

#- get time
timeNow="$(date +"%Y-%m-%d-%H:%M")"

#- create log direct
logDir="compile_log"
cd $MSP_CFDEM_ETC_DIR
mkdir -p $logDir

#- define variables
#- 这里 BASH_SOURCE[0] 等价于 BASH_SOURCE， 表示获取当前执行的 shell 文件所在的路径及文件名
logPath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"
logFileName="log_compile_LIGGGHTS"
logFileContentHeader="$logFileName""_""$timeNow"
doClean="true"

#================================================================================#
# compile LIGGGHTS libraries (forces clean, and then compile)
#================================================================================#
bash $MSP_CFDEM_ETC_DIR/compile_LIGGGHTS_lib.sh
if [ ${PIPESTATUS[0]} -ne 0 ]; then
  exit 1
fi

#================================================================================#
# compile LIGGGHTS src
#================================================================================#
compileLIGGGHTS $logPath $logFileName $logFileContentHeader $doClean
if [ ${PIPESTATUS[0]} -ne 0 ]; then
  exit 1
fi
