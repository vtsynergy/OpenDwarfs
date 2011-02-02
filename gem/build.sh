#!/bin/bash -x

rm -f *.o
rm -f gemnoui
rm -f phi.out

#gcc -O3 -c -I/Developer/GPU\ Computing/OpenCL/common/inc -Iinclude -DNO_UI *.cl
gcc -O3 -c -I/usr/include/CL -Iinclude -DNO_UI -Wall *.c 
g++ -O3 -c -I/usr/include/CL -Iinclude -DNO_UI -Wall *.cpp

g++ -o gemnoui *.o  -lOpenCL
