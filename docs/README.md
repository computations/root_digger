# Root Digger

Root digger is a program that will, given an MSA, tree, and non-reversable
model, place a root on an unrooted tree.

# Building

Root digger depends on `pll-modules`, which can be satisfied by

```
git submodule update --init --recursive
```

Furthermore, the particular variant of `libpll` that is used by Root digger
requires The GNU Scientific Library (GSL). Most distributions have packages
available, so installing it should be done through your package manager.

Root digger requires `cmake` to build. There is a `makefile` provieded which
will set up the build directory and build the software automatically. Once the
software is built, the binary `rd` is placed in the `bin` directory, along with
`rd_test`, which is the test suite.

# Usage

```
./rd --msa <MSA FILE> --tree <TREE FILE> --model <MODEL FILE>
```

The MSA file can be in any format that is supported by `libpll`, which at the
time of writing is one of: "relaxed" phylip or fasta. The tree file should
contain a metric tree in newick format. The model file should be a comma
separated list of numbers which specifies the subsitution model. For example, if
you have some DNA data, and your substitution rate matrix looks like this:


```
*, a, b, c,
d, *, e, f,
g, h, *, i,
j, k, l, *,
```


then the file should have the contents

```
a, b, c, d, e, f, g, h, i, j, k, l
```
