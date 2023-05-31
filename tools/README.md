# VeryFastTree Tools

**treecmp.py**: compare two phylogenetics trees and calculate the Robinson-Foulds distance between them

    python3 treecmp.py tree1 tree2

**prof**: perform profiling of memory and CPU usage

    ./prof FREQ FILE CMD [ARGS]

    FREQ: Time in seconds between each measurement. It can be a decimal value
    FILE: Represent the file where the output will be written
    CMD: The command to be executed
    ARGS: Additional arguments for the command

**vftsum.py**: summarize the stages in the FastTree or VeryFasTree log

    python3 vftsum.py <file|folder> [-s]

    <file|folder>: Log file or folder containing multiple logs with the .out extension
    -s: Perform a shorter summary

    format: VIRT | RES | SHR | MEN% | CPU%

**vftplot.py**: generate plots of memory consumption(RES) and CPU usage using from prof output file

    python3 vftsum.py <file|folder> [step]

    <file|folder>: Log file or folder containing multiple profiles with the .prof extension
    step: Specifies the level of point grouping in a single step for plotting purposes

**winbuild.py**: generates Windows release files

    python3 winbuild.py version

    version: Specifie the name to be used in the generated release files

