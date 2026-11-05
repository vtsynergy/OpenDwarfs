#!/bin/bash

# set PATH and LD_LIBRARY_PATH for CUDA/OpenCL installation (may need to be adjusted depending on where CUDA/OpenCL is installed)
export PATH=$PATH:/usr/local/cuda/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib:/usr/local/cuda/lib64

# call script to set environment variables for HMPP (may need to be adjusted depending on where HMPP is installed)
source ~/HMPPWorkbench-2.5.1/bin/hmpp-env.sh

for currDir1 in *
do
    if [ -d $currDir1 ]
    then
	cd $currDir1
	pwd
	for currDir2 in *
	do
	    if [ -d $currDir2 ]
	    then
		cd $currDir2
		echo $currDir2
		for currDir3 in *
		do
		    if [ -d $currDir3 ]
		    then
			cd $currDir3
			pwd
			make clean
			make exe
			make install
			cd ..
		    fi
		done
		cd ..
	    fi
	done
	cd ..
    fi
done
