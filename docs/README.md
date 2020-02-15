# Root Digger

RootDigger is a program that will, when given a MSA and an unrooted tree place a
root on the given tree.

# Building

Root digger depends on `pll-modules`, which can be satisfied by

    git submodule update --init --recursive

Furthermore, the particular variant of `libpll` that is used by Root digger
requires the GNU Scientific Library (GSL). Most distributions have packages
available, so installing it should be done through your package manager. In the
case that GSL is not found though, the build process will automatically download
GSL and build it.

Compilation requires the dependencies Flex and Bison, both of which can also be
obtained through a package manager.

Root digger requires `cmake` to build. There is a `makefile` provided which will
set up the build directory and build the software automatically. Once the
software is built, the binary `rd` is placed in the `bin` directory, along with
`rd_test`, which is the test suite.

# Usage

    ./rd --msa <MSA FILE> --tree <TREE FILE>

The MSA file can be in any format that is supported by `libpll`, which at the
time of writing is one of: "relaxed" phylip or fasta. The tree file should
contain a metric tree in newick format.

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
