# VeryFastTree

**VeryFastTree** is a highly-tuned implementation of the [FastTree-2](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490) tool that takes advantage of parallelization and vectorization strategies to speed up the inference of phylogenies for huge alignments. It is important to highlight that **VeryFastTree** keeps unchanged the phases, methods and heuristics used by FastTree-2 to estimate the phylogenetic tree. In this way, it produces trees with the same topological accuracy than FastTree-2. In addition, unlike the parallel version of FastTree-2, VeryFastTree is deterministic.

Regarding the performance, for example, **VeryFastTree** (v3.0 - May 2020) is able to construct a tree on a standard server (12-core Intel Xeon E5-2680v3 processor and 128 GiB of memory) using double precision arithmetic from an [ultra-large 330k alignment](http://www.microbesonline.org/fasttree/) in only 4.5 hours, which is 7.8× and 3.5× faster than the sequential and best parallel FastTree-2 times, respectively.

To facilitate the adoption from the research community, VeryFastTree keeps exactly the same command line arguments than FastTree-2. In this way, it is only necessary to replace the call to FastTree-2 by a call to VeryFastTree using the same options to increase the overall performance.

**Release Notes**:

- v3.3.0 (in development):
	- Improved performance with new parallel regions (Local bootstraps, ML Lengths, LogLk, ML splits, etc.).
	- Tree creation (Top hits, TopHitNJSearch, and FastNJSearch, ExhaustiveNJSearch(-slow)) now use threads.
	- Improved non-deterministic mode and no longer uses mutex. 
	- Deterministic mode now also parallelizes non-deterministic parts, but require more computation.
	- New thread levels have been introduced.
	- Tree partitioning method logging now is hidden by default.
	- Tree partitioning algorithm is faster and has a cache.
	- Added FastTree-2.11 changes.
	- Fastq format and libBZ2 compression is now supported. 

- v3.2.0 (December 2022):
    - Decrease in the peak memory usage.
    - Now profiles can be optionally stored on disk. It causes an important reduction in the memory usage.
    - All supported input files (NEXUS, Fasta and Phylip) supports Zlib compression.

- v3.1.0 (December 2021):
    - NEXUS format is now supported.
    - Parallelization of SPR (Subtree-Prune-Regraft) moves.
    - New tree partitioning algorithm for NNI and SPR.
    - Choose between deterministic and non-deterministic parallelization (best performance).
    - Minor fixes.

- v3.0 (May 2020):
    - AVX2 and AVX512 support.
    - NNI (Nearest-Neighbor-Interchange) parallelization.
    - Parallel computation of posterior distributions for each internal node.
    - Deterministic result.

**VeryFastTree** is now included as package in [Bioconda](https://anaconda.org/bioconda/veryfasttree).

If you use **VeryFastTree**, please cite:

[VeryFastTree: speeding up the estimation of phylogenies for large alignments through parallelization and vectorization strategies](https://doi.org/10.1093/bioinformatics/btaa582)  
César Piñeiro, José M. Abuín and Juan C. Pichel.  
Bioinformatics, vol. 36, no. 17, pages 4658-4659, 2020.

# Getting started #

## Requirements

All the libraries needed to compile **VeryFastTree** are included inside the lib folder. The
other basic requirements are:

* CMake v3.5+
* C++11 compiler
    * GCC 5+ (GCC 4 is [bugged](https://gcc.gnu.org/bugzilla/show_bug.cgi?id=56859))
    * Visual studio 2015 (previous versions with support for C++11 may work)
* make (linux only)

## Configuring

CMake will generate the necessary configuration files for the compilation. In Linux/MacOS by
default a native compilation is done for the considered architecture. This
allows to detect automatically the hardware features. If the
program is going to be executed on a different machine, this option must be disabled.

Options can be listed with cmake:

    cmake -LH

    // enable/disable system's processor architecture optimization (linux)
    USE_NATIVE:BOOL=ON

    // enable/disable SSE2 in Windows (Linux default)
    USE_SEE2:BOOL=ON

    // enable/disable SSE4.1
    USE_SEE4:BOOL=OFF

    // enable/disable AVX
    USE_AVX:BOOL=OFF

    // enable/disable AVX2
    USE_AVX2:BOOL=OFF

    // enable/disable AVX512
    USE_AVX512:BOOL=OFF

Example:

    cmake -DUSE_NATIVE=OFF -DUSE_SEE4=ON . // Disable native compilation and use SSE 4.1

## Building

Windows:

    cmake [options] .
    cmake --build

Linux/MacOS:

    cmake [options] .
    make
    make install #Optional


## Running VeryFastTree ##

**To improve the usability and facilitate the adoption of VeryFastTree, it implements the same command interface than FastTree-2. It means that arguments have exactly the same behavior as in FastTree-2.** All these arguments can be consulted with the *-h* option. As a consequence, to take advantage of the performance benefits provided by VeryFastTree is only necessary to replace the call to FastTree-2 by a call to VeryFastTree using the same options.  

VeryFastTree accepts alignments in NEXUS, Fasta, Fastq or Phylip interleaved formats compressed with ZLib and libBZ2.

On the other hand, VeryFastTree has its own extra arguments which have been grouped in the  *Optimizations* section. These arguments are related to the parametrization of the different parallelization and vectorization strategies included in VeryFastTree. Next we list and explain the new arguments available:

- **-threads [n]**
Number of threads (*n*) used in the parallel execution. If this option is not set, the corresponding value will be obtained from the environment variable *OMP\_NUM\_THREADS*. This is the same approach followed by FastTree-2. If *n=1*, VeryFastTree behaves in the same way than FastTree-2 compiled without the *-DOPENMP* flag.

- **-threads-level [level]**
Degree of parallelization: 
	- If level is *0*, VeryFastTree uses the same parallelization strategy as FastTree-2 with some new parallel blocks. 
	- If level is *1*, VeryFastTree uses its tree partitioning method to perform tree computation in sequential order. 
	- If level is *2*, VeryFastTree accelerates the rounds of ML NNIs using its tree partitioning method. 
	- If level is *3* (default), VeryFastTree performs the computation without preserving sequential order.
	- If level is *4*, VeryFastTree accelerates the rounds of SPR steps using its tree partitioning method (it can only be used with datasets larger than 2^sprlength). 

    Note: Each level includes the previous ones, and computation at level *2* and above is performed in a different tree traverse order, so the result may change but is still correct.

- **-threads-mode [mode]**
Changes the mode of parallelization: 
	- If level is *0*, VeryFastTree uses all parallel parts of FastTree-2 including non-deterministic. 
	- If level is *1* (default), VeryFastTree only uses deterministic parallelization parts.

- **-thread-subtrees [num\_subtrees]**
It sets a maximum number of subtrees assigned to the threads. This option could increase the accuracy for small datasets containing large sequences at the expense of reducing the workload balance among threads.

- **-threads-verbose**
To show subtrees assigned to the threads and theoretical speedup, only with verbose > 0

- **-double-precision**
To use double precision arithmetic. Therefore, it is equivalent to compile FastTree-2 with *-DUSE\_DOUBLE*.

- **-ext [type]**
It enables the vector extensions:
    - **AUTO**: (default) selects AVX2 when *-double-precision* is used and SEE3 otherwise. If one extension is not available, the previous level is used.
    - **NONE**: Operations are performed with the native programming language operators. In addition, loops are unrolled with the aim of providing hints to the compiler for applying some optimization (including vectorization).
    - **SSE3**: Arithmetic operations are performed using SSE3 vector intrinsics. Each instruction operates on 128 bit registers, which could contain four 32-bit floats or two 64-bit doubles.
    - **AVX**: Arithmetic operations are performed using AVX vector intrinsics. Each instruction operates on 256 bit registers, which could contain eight 32-bit floats or four 64-bits doubles.
    - **AVX2**: Similar to AVX, but some arithmetic operations are performed using  additional AVX2 vector intrinsics not included in the AVX instruction set. Each instruction operates on 256 bit registers, which could contain eight 32-bit floats or four 64-bit doubles).
    - **AVX512**: Arithmetic operations are performed using AVX512 vector intrinsics. Each instruction operates on 512 bit registers, which could contain sixteen 32-bit floats or eight 64-bits doubles.

- **-fastexp [implementation]**  
This option is used to select an alternative implementation for the exponential function (e<sup>x</sup>), which has a significant impact on performance:
    - **0**: (default) Use the *exp* function included in the built-in math library with double precision.
    - **1**: Use the *exp* function included in the built-in math library with simple precision (not recommended together with *-double-precision* option).
    - **2**: Use a very efficient and fast implementation to compute an accurate approximation of e<sup>x</sup> using double precision arithmetic.
    - **3**: Use a very efficient and fast implementation to compute an accurate approximation of e<sup>x</sup> using simple precision arithmetic (not recommended together with *-double-precision* option).

- **-disk-profiles [0 - 1]**
Reduce the amount of memory required using the disk. Some sequence profiles are stored on disk instead of in memory. Ratio between 0 and 1 specifies the number of profiles to be stored on disk. That is, 1 will store all profiles on disk. Storing more profiles on disk increases the running time. Note that even if the disk file is created, profiles will only be saved to disk when the memory is nearly full.

- **-disk-profiles-file [path]**
Select an alternative file location for storing profiles with **-disk-profiles**. The disk must have enough space to store the profiles. By default the path is the working directory.

- **-disk-profiles-opt**
This parameter forces that all attributes within the profiles are saved to disk. It saves more memory but the file disk could be huge.
