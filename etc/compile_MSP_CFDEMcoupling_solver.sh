#!/bin/bash

#===================================================================#
# Description: Compile MSP-CFDEMcoupling solvers
# Author: Wang Ding
#===================================================================#

source $MSP_CFDEM_TOOLS_DIR/compile_msp_cfdem_solver.sh

# create log direct
logDir="compile_log"
cd $MSP_CFDEM_ETC_DIR
mkdir -p $logDir

#================================================================================#
# compile solver
#================================================================================#
compiledSolverList="$MSP_CFDEM_ETC_DIR/package_msp_cfdem_solver_list.txt"
echo "Please provide the solvers to be compiled in the $compiledSolverList file."

if [[ ! -f $compiledSolverList ]]; then
  echo "$compiledSolverList does not exist. Nothing will be done."
else
  # loop compiledSolverList
  # 修改内部字段分隔符
  oldIFS=$IFS
  IFS=$'\n'
  for line in $(cat $compiledSolverList); do
    if [[ "$line" = solver_dir:* ]]; then
      # 从每一行的信息中获取到要编译的solver所在的目录名
      lineLen=$(expr length "$line")
      solverLen=$(expr $lineLen - 11)
      solverName=$(expr substr "$line" 12 $solverLen)
      path="$MSP_CFDEM_SOLVERS_DIR/$solverName"
      if [ -d $path ]; then
        # compile solver
        cd $path
        logPath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"
        logFileName="log_compile_msp_$solverName"
        logFileContentHeader="$logFileName-$(date +%Y-%m-%d-%H:%M)"
        doClean=true
        compileMspCfdemSolver $path $logPath $logFileName $logFileContentHeader $doClean
        if [ ${PIPESTATUS[0]} -ne 0 ]; then
          IFS=$oldIFS
          exit 1
        fi
      fi
    fi
  done
  IFS=$oldIFS
fi
