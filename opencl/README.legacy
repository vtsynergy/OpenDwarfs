OpenDwarfs
==========

  The OpenDwarfs project provides a benchmark suite consisting of different computation/communication idioms, i.e., dwarfs, for state-of-art multicore CPUs, GPUs, Intel MICs and Altera FPGAs. 
  
  The first instantiation of the OpenDwarfs has been realized in OpenCL, as briefly described in "OpenCL and the 13 Dwarfs: A Work in Progress" by Feng, Lin, Scogland, and Zhang in the 3rd ACM/SPEC International Conference on Performance Engineering, April 2012. 
  The current version, which contains an in-depth performance evaluation on a subset of OpenDwarfs, is described in "On the Characterization of OpenCL Dwarfs on Fixed and Reconfigurable Platforms" by Krommydas, Feng, Owaida, Antonopoulos, and Bellas in the 25th IEEE International Conference on Application-specific Systems, Architectures and Processors (ASAP), June 2014.
  A more thorough description of the latest version, including further in-depth performance evaluation for a larger number of OpenDwarfs, is described in "OpenDwarfs: Characterization of Dwarf-based Benchmarks on Fixed and Reconfigurable Architectures" by Krommydas, Feng, Antonopoulos, and Bellas in Journal of Signal Processing Systems (JSPS), Springer, October 2015.
  
The computation/communication idioms are based on the 13 Berkeley Dwarfs:
(http://view.eecs.berkeley.edu/wiki/Dwarf_Mine).

Benchmark status
----------------

Stable:
gem

Beta:
bfs
cfd
crc
fft
kmeans
lud
nw
spmv
srad
swat
bwa_hmm
nqueens

Alpha:
tdm

Requirements
------------

Packages and libraries needed to build and run the applications.

To build:

    opencl >= 1.0 (some apps require 1.1, but we do not yet guarantee support for 1.2 in all applications.)
    autoconf >= 2.63
    autoheader
    automake
    libtool
    gcc
    maker

To run:

    opencl libs

Building
--------

To build all of the included applications:

    $ ./autogen.sh
    $ mkdir build
    $ cd build
    $ ../configure
    $ make

To build only the applications you select, call configure with the --with-apps
option:

    $ ../configure --with-apps=srad,gem,cfd

To see a full list of options and applications:

    $ ../configure --help

Running
-------

See the application-specific README file in each application's directory.
All the dwarf applications support a common list of options for optionally specifying the OpenCL platform ID (-p)
and OpenCL device ID (-d), or alternatively, the device type (-t). Optionally you can provide -o option to use 
optimized kernels. It picks up the optimized kernel for the given device type. For an example, if the device in use 
is GPU and -o option is provided, it will use <kenel_name>_opt_gpu.cl file present in the application directory. 
These options, if supplied, must follow the executable name and be delimited from the application-specific options by double dashes (--).

General format: ./<executable> [-p <platform> -d <device> | -t <type> -o --] [app-specific options]

    <platform>	:integer ID of platform to use
    <device>    :integer ID of device in <platform> to use
    <type>	    : device type to use (0:CPU, 1:GPU, 2:MIC, 3:FPGA)
    -o          :Optional flag to use the optimzed flag for the device in use

Example1: ./astar -p 0 -d 0 -- (selects device with device ID 0 on platform with platform ID 0)
Example2: ./astar -t 0 -- (selects CPU device type on default platform with platform ID 0, if available)
Example2: ./nw -p 0 -d 0 -o -- (Run the optimized dwarf of device ID 0 on platfrom with platform ID 0)

Notes:	If no parameters are supplied, default platform ID is 0 and default device type is CPU.
	If -t parameter is given, default platform ID 0 is searched for supplied device type <type>. If not available, CPU device type selection will be attempted.
	If device ID is unknown, a combination of -p and -t is available to search for device of selected <type> on platform ID <platform>.
    If the optimized kernel does not exist, application wil throw and error and exit. 

Notes: SWAT DOES NOT compile for OpenCl and FFT kernel DID NOT fit on Stratix V in this release. 

Acknowledgements
----------------

This project has been supported in part by Air Force Research Lab, Altera, AMD, Department of Defense, Harris, Los Alamos National Laboratory, and Xilinx via the NSF Center for High-Performance Reconfigurable Computing (CHREC) under NSF grant IIP-0804155 and indirectly by AFOSR grant FA9550-12-1-0442 and NSF grants CNS-0916719 and MRI-0960081.

Integration for Altera FPGA support for crc and csr, as well as extensions for these
benchmarks, have been contributed by Tyler Kenney at IBM.

Part of the OpenDwarfs benchmark suite (as acknowledged in the respective benchmarks' READMEs) was ported to OpenCL from the corresponding CUDA implementations in earlier implementations of the Rodinia benchmark suite (http://www.cs.virginia.edu/~skadron/wiki/rodinia/index.php/Main_Page). 
