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

VeryFastTree preserves all legacy options from FastTree, new options perform optimizations 
that greatly increase performance.
