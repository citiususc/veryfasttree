# VeryFastTree

# Getting started #

## Requirements

All libraries need to compile **VeryFastTree** are include inside of lib folder, the 
other basic requirement are:

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

VeryFastTree implements through the CLI library the same command interface that FastTree used. The arguments have the same behavior as in FastTree and these can be consulted with *-h* option . VeryFastTree adds new arguments grouped in optimizations section, the arguments allow to improve the performance and to change between the different versions that in FastTree were different binaries.

  - **-threads**: allows to specify the number of threads, if the option
    is not set, value will be obtained from the environment variable
    *OMP\_NUM\_THREADS* as FastTree. If the number of threads is 1, the
    tool behavior is the same as FastTree compiled without -DOPENMP.

  - **-threads-balanced**: Removes the maximum number of subtrees
    restriction assigned to each thread, improves the performance of the
    parallel version but reduces the interaction between nodes when
    creating more isolated zones, it can reduce the accuracy with some
    datasets.

  - **-threads-level**: allows to change the type of parallelization, 0 is
    the parallelization used by FastTree and 1 (by default) the new
    parallelization system where tree processing is divided in threads.

  - **-double-precision**: Use double precision numbers for operations.
    The same behavior as compiling FastTree with *-DUSE\_DOUBLE*.

  - **-fastexp**: allows selecting an alternative implementation of the
    exp function, which has a significant impact on performance:
    
      - 0: (default) use math exp with double precision.
    
      - 1: use math exp with simple precision. (not recommended with
        -double-precision option)
    
      - 2: Approximate vectorized exp from Cephes Math Library with
        double precision..
    
      - 3: Approximate vectorized exp from Cephes Math Library with
        simple precision.. (not recommended with -double-precision
        option)

  - **-ext**: Enable vector extensions:
    
      - none: (default) Operations are performed with language
        operators, loops are unrolled which allows the compiler to apply
        some optimization, including some type of vectorization.
    
      - SSE, AVX128: Use the set of 128-bit vector instructions, 4
        floats of 32 bits or 2 doubles of 64 bits, it requires data
        aligning which a small memory overhead.
    
      - AVX256: Use the set of 256-bit vector instructions, 8 floats of
        32 bits or 4 doubles of 64 bits, doubles require data aligning
        which a small memory overhead and floats add a padding to
        preserve aligning when use 8 that is not divisor of protein set
        which requires more memory.


