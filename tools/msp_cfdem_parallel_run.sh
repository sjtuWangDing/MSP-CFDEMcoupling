#!/bin/bash:

#===================================================================#
# Description: Defining function to run a parallel CFD-DEM case
# Author: Wang Ding
#===================================================================#

mspCfdemParallelRun() {
  # define variables
  casePath="$1"
  solver="$2"
  numberOfProcs="$3"
  machineFileName="$4"
  logPath="$5"
  logFileName="$6"
  decomposeCase="$7"

  # decompose case
  if [[ $decomposeCase == "false" ]]; then
    echo "decomposePar: Not decomposing case."
  else
    echo "decomposePar: Decomposing case."
    cd $casePath/CFD
    decomposePar -force
  fi

  # make processor dirs visible
  for i in $(seq 1 $numberOfProcs); do
    count=$[ $i - 1 ]
    (cd $casePath/CFD/processor$count && touch file.foam)
  done

  # clean old log file
  if [ -f "$logPath/$logFileName" ]; then
    rm $logPath/$logFileName
  fi

  # 重定向，将文件描述符2（标准错误输出）重定向到 1（STDOUT_FILENO）
  # 如果既想把输出保存到文件中，又想在屏幕上看到输出内容，就可以使用 tee 命令
  # tee -a file: 输出到标准输出的同时，追加到文件file中。如果文件不存在，则创建；如果已经存在，就在末尾追加内容，而不是覆盖
  # 如果不加 &1 直接 echo 2>1  就变成重定向输出到 1 这个文件里去了，如果没有 1，系统就自动创建一个文件 1
  echo "\
/*---------------------------------------------------------------------------*\\
|                            Run Log Content                                 |
\*---------------------------------------------------------------------------*/
" 2>&1 | tee -a /$logPath/$logFileName
  echo -n "Log path: " 2>&1 | tee -a /$logPath/$logFileName
  echo "$logPath/$logFileName" 2>&1 | tee -a /$logPath/$logFileName
  echo

  if [ "none" = $machineFileName ]; then
    mpirun -np $numberOfProcs $solver -parallel 2>&1 | tee -a $logPath/$logFileName
  else
    mpirun -machineFile $machineFileName -np $numberOfProcs $solver -parallel 2>&1 | tee -a $logPath/$logFileName
  fi
}
