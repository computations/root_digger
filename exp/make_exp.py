#!/usr/bin/env python3

import argparse
import csv
import datetime
import math
import os
import random
import shutil
import subprocess
import sys

import ete3
import numpy

if not shutil.which("indelible"):
    print("Please add indelible to your path")
    sys.exit()

if not shutil.which("iqtree"):
    print("Please add iqtree to your path")
    sys.exit()

CONTROL_FILE= """
[TYPE] NUCLEOTIDE 1

[SETTINGS]
  [randomseed] {randomseed}

[MODEL] m1
  [submodel] UNREST {model_params}
  [statefreq] {freq_params}

[TREE] t1 {tree}

[PARTITIONS] p1
  [t1 m1 {sites}]

[EVOLVE] p1 1 seqs
"""

RD = os.path.abspath("../bin/rd") + " --msa {msa} --tree {tree} --states 4 --seed {seed} --silent --slow"
model_file = "subst.model"
freqs_file = "freqs.model"

IQTREE_RF = "iqtree -rf_all {trees}"
IQTREE_R  = "iqtree -seed {random_seed} -r {taxa} {outfile}"

TAXA_STEPS = [10, 100]
SITE_STEPS = [100, 1000]
RUN_TEMPLATE = "run_{run_iter:0{leading_zeroes}}"
TOTAL_ITERS = 100

class directory_guard:
    def __init__(self, path):
        self._path = path
    def __enter__(self):
        self._old_dir = os.getcwd()
        os.chdir(self._path)
        return self
    def __exit__(self, *args):
        os.chdir(self._old_dir)

class subst_params:
    def __init__(self):
        self._params = numpy.random.rand(4,4)
        self._params -= numpy.diag(numpy.diag(self._params))
        self._params -= numpy.diagflat(numpy.dot(self._params, 
            numpy.ones((4,1))))
    def indel_repr(self):
        p = []
        for i in range(4):
            for j in range(4):
                if i == j:
                    continue
                p.append(self._params[i][j])
        return ' '.join([str(f) for f in p])
    def rd_repr(self):
        P = numpy.array([[0,0,0,1],
                         [0,1,0,0],
                         [1,0,0,0],
                         [0,0,1,0]])
        tmp = numpy.inner(numpy.inner(P.T,self._params), P)
        p = []
        for i in range(4):
            for j in range(4):
                if i == j:
                    continue
                p.append(tmp[i][j])
        return ','.join([str(f) for f in p])

class freq_params:
    def __init__(self):
        self._params = numpy.random.dirichlet([1.0 for _ in range(4)])
    def indel_repr(self):
        return ' '.join([str(f) for f in self._params])
    def rd_repr(self):
        P = numpy.array([[0,0,0,1],
                         [0,1,0,0],
                         [1,0,0,0],
                         [0,0,1,0]])
        p = numpy.inner(self._params, P)
        return ','.join([str(f) for f in p])

class exp:
    def __init__(self, root_path, run_iter, seed=None):
        leading_zeroes = math.ceil(math.log10(TOTAL_ITERS))
        self._run_path = os.path.abspath(os.path.join(root_path,
            RUN_TEMPLATE.format(run_iter= run_iter, leading_zeroes =
                leading_zeroes)))
        if not os.path.exists(self._run_path):
            os.mkdir(self._run_path)

        self._seed_file = os.path.join(self._run_path, '.seed')
        if not os.path.exists(self._seed_file):
            if seed is None:
                self._seed = int.from_bytes(os.urandom(4), 'little')
            else:
                self._seed = seed
            with open(self._seed_file, 'w') as sf:
                sf.write(str(self._seed))
        else:
            with open(self._seed_file) as sf:
                self._seed = int(sf.read())

    @staticmethod
    def check_done_indel(path):
        if os.path.exists(os.path.join(path, '.done')):
            with open(os.path.join(path, ".done")) as done_file:
                for line in done_file:
                    if line.find('indel') != -1:
                        return True
        return False

    @staticmethod
    def set_indel_done(path):
        with open(os.path.join(path, ".done"), 'a') as done_file:
            done_file.write("indel:" + datetime.datetime.now().isoformat())

    @staticmethod
    def check_done_rd(path):
        if os.path.exists(os.path.join(path, '.done')):
            with open(os.path.join(path, ".done")) as done_file:
                for line in done_file:
                    if line.find('rd') != -1:
                        return True
        return False

    @staticmethod
    def set_rd_done(path):
        with open(os.path.join(path, ".done"), 'a') as done_file:
            done_file.write("rd:" + datetime.datetime.now().isoformat())

    def get_model_params(self):
        numpy.random.seed(self._seed)
        return (freq_params(), subst_params())

    def make_indel_control_file(self, freqs, subst, tree, sites):
        return CONTROL_FILE.format(
                freq_params = freqs.indel_repr(),
                model_params = subst.indel_repr(),
                tree = tree,
                sites = sites,
                randomseed = self._seed)

    def run_exp(self, sites, freqs, subst, tree_filename):
        with open(tree_filename) as treef:
            tree = treef.read()
        with open('control.txt', 'w') as control_txt:
            control_txt.write(self.make_indel_control_file(freqs, subst, tree,
                sites))
        if not self.check_done_indel('.'):
            subprocess.run("indelible", stdout=subprocess.DEVNULL)
            self.set_indel_done('.')

        if not self.check_done_rd('.'):
            rd_output = subprocess.run(RD.format(msa="seqs_TRUE.phy",
                tree=os.path.join("../", tree_filename),
                seed=self._seed).split(' '),
                stdout=subprocess.PIPE)
            with open('rd_output', 'w') as rd_outfile:
                rd_outfile.write(rd_output.stdout.decode('utf-8'))
            self.set_rd_done('.')

    def run_all(self):
        old_dir = os.getcwd()
        os.chdir(self._run_path)

        freqs, subst = self.get_model_params()
        with open('subst.model', 'w') as model_file:
            model_file.write(subst.rd_repr())

        with open('freqs.model', 'w') as model_file:
            model_file.write(freqs.rd_repr())

        for taxa in TAXA_STEPS:
            all_trees = []
            tree_file = os.path.abspath(os.path.join('../..',str(taxa)+".tree"))
            with open(tree_file) as tf:
                all_trees.append(tf.read())
            for sites in SITE_STEPS:
                exp_dir = "{taxa}taxa_{sites}sites".format(taxa=taxa,
                            sites=sites)
                if not os.path.exists(exp_dir):
                    os.mkdir(exp_dir)
                with directory_guard(exp_dir):
                    self.run_exp(sites, freqs, subst, tree_file)
                    with open('rd_output') as tf:
                        all_trees.append(tf.read())
            result_tree_file = "result_trees_{}_taxa".format(taxa)
            with open(result_tree_file, 'w') as rt:
               rt.write(''.join(all_trees))

            parsed_trees = [ete3.Tree(t) for t in all_trees]
            true_tree = parsed_trees[0]

            rfdists = []
            for i in range(1, len(parsed_trees)):
                rfdists.append(compute_root_distance(true_tree,
                    parsed_trees[i]))
            with open('rfdists_{taxa}_taxa'.format(taxa=taxa), 'w')\
                    as rf_outfile:
                rf_outfile.write(','.join([str(i) for i in SITE_STEPS]))
                rf_outfile.write('\n')
                rf_outfile.write(','.join([str(f) for f in rfdists]))
                rf_outfile.write('\n')
        os.chdir(old_dir)

def compute_root_distance(t1, t2):
    result = t1.robinson_foulds(t2)
    return result[0]

def summarize_results(path):
    with directory_guard(path):
        for taxa in TAXA_STEPS:
            leading_zeroes = math.ceil(math.log10(TOTAL_ITERS))
            nrows = len(SITE_STEPS)
            totals = numpy.zeros((nrows,))
            for i in range(TOTAL_ITERS):
                result_tree_file = os.path.join(RUN_TEMPLATE.format(run_iter=i,
                    leading_zeroes = leading_zeroes),
                        "rfdists_{taxa}_taxa".format(taxa=taxa))
                with open(result_tree_file) as rfdist_file:
                    rfdists = csv.DictReader(rfdist_file)
                    for row in rfdists:
                        for i in range(len(SITE_STEPS)):
                            totals[i] = row[str(SITE_STEPS[i])]
            totals /= TOTAL_ITERS
            with open('summary_{}_taxa'.format(taxa), 'w') as outfile:
                outfile.write(','.join([str(i) for i in SITE_STEPS]))
                outfile.write('\n')
                outfile.write(','.join([str(f) for f in totals]))
                outfile.write('\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, help='Path to store the exp',
            required=True)
    args = parser.parse_args()
    exp_path = os.path.abspath(args.path)

    if not os.path.exists(exp_path):
        os.mkdir(exp_path)

    with directory_guard(exp_path):
        for i in range(TOTAL_ITERS):
            print("[", datetime.datetime.now().isoformat(), "]", sep = '', end = '')
            print(" trial:", i)
            e = exp('.', i)
            e.run_all()
        summarize_results('.')
