#!/bin/bash

echo "Vary state with S = 2, T = 1000"
echo -e "State\tTime (s)\tLog_likelihood"
for (( i = 0; i < 4 ; i++ ))
do
    mul=1
    for (( k = 0 ; k < i ; k++ ))
    do
        mul=`expr $mul \* 10`
    done
    for (( j = `expr $mul \* 1`; j < `expr $mul \* 10`; j = `expr $j + $mul` ))
    do
        ./bwa_opencl -n $j -v n
    done
done

echo "Vary symbols with N = 60, T = 1000"
echo -e "Symbols\tTime (s)\tLog_likelihood"
for (( i = 0; i < 4 ; i++ ))
do
    mul=1
    for (( k = 0 ; k < i ; k++ ))
    do
        mul=`expr $mul \* 10`
    done
    for (( j = `expr $mul \* 1`; j < `expr $mul \* 10`; j = `expr $j + $mul` ))
    do
        ./bwa_opencl -s $j -v s
    done
done

echo "Vary observations with N = 60, T = 2"
echo -e "Observations\tTime (s)\tLog_likelihood"
for (( i = 10; i < 100 ; i = `expr $i + 10` ))
do
    ./bwa_opencl -t $i -v t
done
for (( i = 100; i < 1000 ; i = `expr $i + 100` ))
do
    ./bwa_opencl -t $i -v t
done
for (( i = 1000; i < 10000 ; i = `expr $i + 1000` ))
do
    ./bwa_opencl -t $i -v t
done


