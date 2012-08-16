Overview
========

This repository contains the source code for the BITS (Binary Interval Search)
algorithm developed in the Quinlan laboratory at the University
of Virginia.

Installation
============
BITS installation requires the GNU Scientific Library (GSL) as well as
NVidia GPU drivers and libraries. Below is a step-by-step tutorial for how
to install the necessary drivers, SDK and libraries to run BITS on an
NVidia CUDA GPU.  If you have questions, email me.

Install the GNU Scientific Libraries (GSL).

    - This is typically quite simple, as one can use package managers.
        - e.g., for OS X using Homebrew: brew install gsl
        - e.g., for Ubuntu: apt-get gsl

Download and install the CUDA toolkit.

    - The CUDA toolkit for multiple platforms is available at:
        - http://developer.nvidia.com/cuda-downloads
    - e.g., for OS X: http://developer.download.nvidia.com/compute/cuda/4_2/rel/toolkit/cudatoolkit_4.2.9_macos.pkg

Download and install the CUDA drivers for you system.

    - Likewise available for many platforms at: 
        - http://developer.nvidia.com/cuda-downloads
    - e.g., for OS X: http://developer.download.nvidia.com/compute/cuda/4_2/rel/drivers/devdriver_4.2.10_macos.dmg

Download and install the CUDA SDK.

    - Likewise available for many platforms at: 
        - http://developer.nvidia.com/cuda-downloads
    - e.g., for OS X: http://developer.download.nvidia.com/compute/cuda/4_2/rel/sdk/gpucomputingsdk_4.2.9_macos.pkg
    - *Important*: Take note of the installation path - you will need this when we update the BITS Makefile.
    - The default path for OSX is `/Developer/GPU\ Computing/`.

Download and install the cudapp library at: http://code.google.com/p/cudpp/. Once downloaded, do
::

        cd cudpp_src_2.0
        cmake .
        make

        *Important*: Take note of the path to which you placed cudpp, as you will need this when we update the BITS Makefile.

Update your PATH (rec. you save in .bash_profile for permanence) as follows
::

        export PATH=/usr/local/cuda/bin:$PATH
        export DYLD_LIBRARY_PATH=/usr/local/cuda/lib:$DYLD_LIBRARY_PATH

Clone the BITS repository
::

        git clone git://github.com/arq5x/bits.git

Navigate into the bits directory
::

        cd bits

Edit the `defs.cuda` file in accordance with your configuration.
    - Edit the `PLATFORM` environment variable.
        * if OSX,   set PLATFORM=darwin
        * if Linux, set PLATFORM=linux
        * if Windows, sorry this is unsupported.
    - Edit the `SDK_PATH` environment variable.
        * This is the installation path that you should have taken note of
          in step #4.
        * e.g., for OS X, this should be: SDK_PATH=/Developer/GPU\ Computing/
    - Edit the `SDK_PATH` environment variable.
        * This is the path to which you downloaded and compiled cudpp in step
          #5.
    - EDIT the `CUDA_LIB` environment variable.
        * If OS X, this should be: /usr/local/cuda/lib
        * If Linux, this should be: /usr/local/cuda/lib64

At this point, you should be ready to compile BITS
::

        make


Now, you can test both the sequential and CUDA versions of the tools by
running the `bits_tests` scripts. Also, this shell script demonstrates how 
to run each of the BITS tools
::

        sh bits_tests

If all works well, you should see the following
::

        sh bits_tests
        bits_count
        72534
        bits_count_cuda
        72534
        bits_count_per_interval
        72534
        bits_count_per_interval_cuda
        72534
        bits_enumerate
        72534
        bits_enumerate_cuda
        72534
        bits_test
        O:72534 E:1124.853000   sd:33.680585    p:0.000999
        bits_test_cuda
        O:72534 E:1124.081000   sd:36.552024    p:0.000999


Usage
=====

If you want to run the BITS CUDA Monte Carlo simulation tool for a single
pairwise comparison of two BED files, run the following, where -n is the 
number of MC iterations and -g is the name and size of each chromosome::

    bin/bits_test_cuda -a a.bed -b b.bed -g chrom.sizes -n 1000 

If you want to run the BITS CUDA Monte Carlo simulation tool for a _many_
pairwise comparison of multiple BED files, just create a shell script that
loops over every pairwise set of files and calls the program as above. For
example::
    
    for file1 in `cat file_list.txt`
    do
        for file2 in `cat file_list.txt`
        do
            bin/bits_test_cuda -a $file1 -b $file2 -g chrom.sizes -n 1000
        done
    done

