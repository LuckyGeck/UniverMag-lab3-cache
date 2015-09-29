#!/usr/bin/env bash
matrixSizes=( 127 128 129 1023 1024 1025 1040 1041 1050 1100 )
cacheSizes=( 32 256 3072 )
cacheLineSizes=( 8 8 12 )
for matrixSize in "${matrixSizes[@]}"
do
    for i in `seq 3`
    do
        echo '$> ./cache' $matrixSize ${cacheSizes[$i-1]} 64 ${cacheLineSizes[$i-1]}
        ./cache $matrixSize ${cacheSizes[$i-1]} 64 ${cacheLineSizes[$i-1]}
    done
done
