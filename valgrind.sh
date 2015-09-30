#!/usr/bin/env bash


matrixSizes=( 127 128 129 1023 1024 1025 1040 1041 1050 1100 )
for matrixSize in "${matrixSizes[@]}"
do
    # valgrind --tool=cachegrind --D1=32768,8,64 --L2=262144,8,64 --LL=3145728,12,64 ./cache_old $matrixSize
    valgrind --tool=cachegrind ./cache_old $matrixSize
done