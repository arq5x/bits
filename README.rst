Overview
========

This repository contains the source code for the BITS (Binary Interval Search)
algorithm developed by Ryan Layer in the Quinlan laboratory at the University
of Virginia.

Installation
============
Installation requires NVidia GPU drivers and libraries.  Stay tuned for
detailed installation instructions.

CUDA download center: http://developer.nvidia.com/cuda/cuda-downloads

0. Install GSL.
	- e.g., for OS X: brew install gsl
	- e.g., for Ubuntu: apt-get gsl

1. Download and install the CUDA toolkit.
    - http://developer.nvidia.com/cuda-downloads
    - e.g., for OS X: http://developer.download.nvidia.com/compute/cuda/4_2/rel/toolkit/cudatoolkit_4.2.9_macos.pkg

2. Download and install the CUDA drivers for you system.
    - http://developer.nvidia.com/cuda-downloads
    - e.g., for OS X: http://developer.download.nvidia.com/compute/cuda/4_2/rel/drivers/devdriver_4.2.10_macos.dmg

3. Download and install the CUDA SDK.
	- http://developer.download.nvidia.com/compute/cuda/4_2/rel/sdk/gpucomputingsdk_4.2.9_macos.pkg
	- default for OSX is /Developer/GPU Computing/

4. Download and install the cudapp library
http://code.google.com/p/cudpp/
  cd cudpp_src_2.0
  cmake .
  make

4. Update PATH (and possibly add to .bash_profile)
	export PATH=/usr/local/cuda/bin:$PATH
	export DYLD_LIBRARY_PATH=/usr/local/cuda/lib:$DYLD_LIBRARY_PATH

4. Edit defs.cuda accordingly.
	a. if OSX,   set PLATFORM=darwin
	   if Linux, set PLATFORM=linux
	   if Windows, sorry.
	
	b. edit SDK_PATH
		- e.g., for OS X, this should be: SDK_PATH=/Developer/GPU\ Computing/
		
	c. edit CUDPP_PATH
		- set this to the path in which you downloaded and compiled cudpp

