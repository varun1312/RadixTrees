Efficient Non-blokcing Radix Trees
----------------------------------

This project contains the source code and proof for the work in the paper titled "Efficient Non-blocking Radix Trees".
This work is accepted for publication at EuroPar 2017 conference (http://europar2017.usc.es/).

Directory Structure
-------------------

The directory structure is organized as follows:

* [code] (./code)
  
The "code" directory contains the source code and the corresponding Makefile.

Build system
------------

This software is tested on a system running on ubuntu 16.04 and RHEL 6.7 systems.

The code is compiled using gcc 4.8.4. Any later version of gcc should work.

Compilation
-----------

The source code has a Makefile corresponding to it. Do "make" in the [code] (./code) directory and an output executable "nbrt.o" is generated

Execution
---------

The executable "nbrt.o" takes 6 command line arguments. These are

***./nbrt.o numOfThreads searchPercentage insertPercentage removePercentage durationInSeconds maximumKeySize initialSeed***

Example: ./nbrt.o 32 70 20 10 50 100000 0

Note that the sum of "searchPercentage", "insertPercentage" and "removePercentage" should amount to 100


License
-------

This software is available under the MIT license.
