# VeryFastTree

# Getting started #

## Requirements

All libraries needed to compile **VeryFastTree** are included inside the lib folder. The
other basic requirements are:

* CMake v3.5+
* C++11 compiler
	* GCC 5+ (GCC 4 is [bugged](https://gcc.gnu.org/bugzilla/show_bug.cgi?id=56859))
	* Visual studio 2015 (previous versions with support for C++11 may work)
* make (linux only)

## Configuring

CMake will generate the necessary configuration files for the compilation. In Linux by
default a native compilation is done for the architecture where it is compiled, this
allows to detect the characteristics automatically and apply them when possible, if the
program is going to be executed in a different machine this option must be disabled.

Options can be listed with cmake:

	cmake -LH

	// enable/disable system's processor architecture optimization (linux)
	USE_NATIVE:BOOL=ON

	// enable/disable SSE2 in windows (linux default)
	USE_SEE2:BOOL=ON

	// enable/disable SSE4.1
	USE_SEE4:BOOL=OFF

	// enable/disable AVX
	USE_AVX:BOOL=OFF

	// enable/disable AVX2
	USE_AVX2:BOOL=OFF

Example:

	cmake -DUSE_NATIVE=OFF -DUSE_SEE4=ON . //Disable native compilation an use SSE 4.1

## Build

Windows:

	cmake [options] .
	cmake --build

Linux:

	cmake [options] .
	make
	make install #Optional


## Running VeryFastTree ##

**To improve the usability and facilitate the adoption of VeryFastTree, it implements the same command interface than FastTree-2. It means that arguments have exactly the same behavior as in FastTree-2.** All these arguments can be consulted with the *-h* option. As a consequence, to take advantage of the performance benefits provided by VeryFastTree is only necessary to replace the call to FastTree-2 by a call to VeryFastTree using the same options.  

On the other hand, VeryFastTree has its own extra arguments which have been grouped in the  *Optimizations* section. These arguments are related to the parametrization of the different parallelization and vectorization strategies included in VeryFastTree. Next we list and explain the new arguments available:
- **-threads [n]**
It allows to specify the number of threads (*n*) used in the parallel execution. If this option is not set, the corresponding value will be obtained from the environment variable *OMP\_NUM\_THREADS*. This is the same approach followed by FastTree-2. If *n=1*, VeryFastTree behaves in the same way than FastTree-2 compiled without the *-DOPENMP* flag.

- **-threads-level [level]**
It allows to change the degree of parallelization. If level is *0*, VeryFastTree uses the same parallelization strategy followed by FastTree-2. If level is *1* (by default), VeryFastTree accelerates the rounds of ML NNIs using its tree partitioning method.

- **-thread-subtrees [num\_subtrees]**
It sets a maximum number of subtrees assigned to a thread. This option could increase the accuracy for small datasets containing large sequences at the expense of reducing the workload balance among threads.

- **-double-precision**
Use double precision arithmetic. Therefore, it is equivalent to compile FastTree-2 with *-DUSE\_DOUBLE*.

- **-ext [type]**
It enables the vector extensions:
	- **none**: (default) Operations are performed with the native programming language operators. In addition, loops are unrolled with the aim of providing hints to the compiler for applying some optimization (including vectorization).
	- **SSE3**: Arithmetic operations are performed using SSE3 vector intrinsics. Each instruction operates on 128 bit registers, which could contain four 32-bit floats or two 64-bit doubles.
	- **AVX**: Arithmetic operations are performed using AVX vector intrinsics. Each instruction operates on 256 bit registers, which could contain eight 32-bit floats or four 64-bits doubles.
	- **AVX2**: Similar to AVX, but some arithmetic operations are performed using  additional AVX2 vector intrinsics not included in the AVX instruction set. Each instruction operates on 256 bit registers, which could contain eight 32-bit floats or four 64-bit doubles).

- **-fastexp [implementation]**  
This option is used to select an alternative implementation for the exponential function (e<sup>x</sup>), which has a significant impact on performance:
	- **0**: (default) Use the *exp* function included in the built-in math library with double precision.
	- **1**: Use the *exp* function included in the built-in math library with simple precision (not recommended together with *-double-precision* option).
	- **2**: Use a very efficient and fast implementation to compute an accurate approximation of e<sup>x</sup> using double precision arithmetic.
	- **3**: Use a very efficient and fast implementation to compute an accurate approximation of e<sup>x</sup> using simple precision arithmetic (not recommended together with *-double-precision* option).
