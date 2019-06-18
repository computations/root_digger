#!/usr/bin/env python3

import os
import shutil
import subprocess
import sys

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

RD = os.path.abspath("../bin/rd") + " --msa {msa} --tree {tree} --model {model} --freqs {freqs}"
model_file = "random_unsym.model"
freqs_file = 'uniform.freqs'

if not os.path.exists('active_exp'):
    os.mkdir('active_exp')

for taxa in [10, 100, 1000, 10000]:
    for sites in [100, 1000, 10000]:
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
        print(RD.format(msa="seqs_TRUE.phy", tree=os.path.join("../..",
            tree_file), model=os.path.join('../../',model_file),
            freqs=os.path.join("../../", freqs_file)).split(' '))
        rd_output = subprocess.run(RD.format(msa="seqs_TRUE.phy", tree=os.path.join("../..",
            tree_file), model=os.path.join('../../',model_file),
            freqs=os.path.join("../../", freqs_file)).split(' '),
            capture_output=True)
        with open('rd_output', 'w') as rd_outfile:
            rd_outfile.write(rd_output.stdout.decode('utf-8'))
        os.chdir(old_dir)
