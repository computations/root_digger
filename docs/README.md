# Root Digger

RootDigger is a program that will, when given a MSA and an unrooted tree with
branch lengths place a root on the given tree. For the foreseeable future,
RootDigger will only support DNA data, as the method RootDigger uses is
ineffective when using AA data.

# Building

Currently, the best way to get the most recent version of RootDigger is by
cloning the repository

    git clone --recursive https://github.com/computations/root_digger

This will obtain all the required dependencies including a modified `libpll`.
Furthermore, the particular variant of `libpll` that is used by RootDigger
requires the GNU Scientific Library (GSL). Most distributions have packages
available, so installing it should be done through your package manager. In the
case that GSL is not found though, the build process will automatically download
GSL and build it.

Compilation of `libpll` requires the dependencies Flex and Bison, both of which
can also be obtained through a package manager.

Root digger requires `cmake` to build. There is a `makefile` provided which will
set up the build directory and build the software automatically. Once the
software is built, the binary `rd` is placed in the `bin` directory, along with
`rd_test`, which is the test suite.

## Parallelism

As of version 1.4, root digger supports a rudimentary version of both thread and
process level parallelism. Both are optional. By default, if cmake can find
OpenMP, it will build with it. 

For MPI, a special build flag needs to be passed. There is a `mpi` target for
the makefile, so  in most cases running `make mpi` should be sufficient. In
cases where this doesn't work, the base command is

    cmake -DMPI_BUILD=ON -DCMAKE_BUILD_TYPE=Release

In general, MPI is better to use on larger trees.

# Usage

    ./rd --msa <MSA FILE> --tree <TREE FILE>

The MSA file can be in any format that is supported by `libpll`, which at the
time of writing is one of: "relaxed" phylip or fasta. The tree file should
contain a metric tree (this is to say, branch lengths which are in expected
substitutions per site) in newick format.

By default, RootDigger runs in search mode with early stopping on. This means
that RootDigger will simply look for the most likely root, and will "stop
early". This means that RootDigger will consider the search concluded when it
finds the same root placement twice in a row (as opposed to requiring that the
likelihood is the same twice in a row).

RootDigger can be run in exhaustive mode using `--exhaustive`. This will cause
RootDigger to consider every branch, and report the Likelihood Weight Ratio of
placing the root on that branch. Informally, this can be interpreted as the
probability of placing the root on a given branch. By default, `--exhaustive`
mode does _not_ run with early stopping. This can be enabled using
`--early-stop`. In practice, this doesn't affect the results at all, but in
principle it could, so be warned.

For more information about the options, there is a `--help` flag which will
print detailed information about all the options.

# TL;DR

**Dependencies**: 

- GSL
- Cmake 3.0
- C++14 compatible compiler

**Optional Dependencies**:

- OpenMP for thread level parallelism
- Some MPI package for process level parallelism

**Usage**:

For search mode:

    rd --msa <MSA> --tree <TREE> 

For exhaustive mode:

    rd --msa <MSA> --tree <TREE> --exhaustive
