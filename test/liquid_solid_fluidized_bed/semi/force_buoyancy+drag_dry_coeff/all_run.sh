#!/bin/bash

currentPath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
for fileName in "$currentPath"/*; do
  if [ -d "$fileName" ]; then
    if [[ "$fileName" != */Yang && "$fileName" != */Dallavalle && "$fileName" != */DiFelice ]]; then
      cd "$fileName"
      ./run.sh > /dev/null
    fi
  fi
done
