#!/usr/bin/env python3

import os
import sys
import shutil
import subprocess

if not shutil.which("indelible"):
    print("Please add indelible to your path")
    sys.exit()

CONTROL_FILE= """
[TYPE] NUCLEOTIDE 1

[MODEL] m1
  [submodel] UNREST 0.541224 0.204832 0.758911 0.164193 0.884369 0.740007 0.405492 0.078174 0.059229 0.420233 0.340426 0.247193

[TREE] t1 {tree}

[PARTITIONS] p1
  [t1 m1 {sites}]

[EVOLVE] p1 1 seqs
"""

if not os.path.exists('active_exp'):
    os.mkdir('active_exp')

for taxa in [10, 100, 1000, 10000]:
    for sites in [100, 1000, 10000, 100000]:
        print("taxa: ", taxa, "sites: ", sites)
        exp_dir = os.path.join("active_exp",
                "{taxa}taxa_{sites}sites".format(taxa=taxa, sites=sites))
        if not os.path.exists(exp_dir):
            os.mkdir(exp_dir)
        tree_file = str(taxa) + ".tree"
        with open(tree_file) as treef:
            tree = treef.read()
        control_file = os.path.join(exp_dir, "control.txt")
        with open(control_file, 'w') as cf:
            cf.write(CONTROL_FILE.format(tree=tree, sites=sites))
        old_dir = os.getcwd()
        os.chdir(exp_dir)
        subprocess.run("indelible", shell=True, check=True)
        os.chdir(old_dir)

