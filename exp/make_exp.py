#!/usr/bin/env python3

import argparse
import csv
import datetime
import math
import multiprocessing
import os
import random
import shutil
import subprocess
import sys
import string

import ete3
import numpy
import progressbar
from Bio import SeqIO
import Bio


progressbar.streams.flush()

PROGRESS_BAR = None
PROGRESS_BAR_ITER = multiprocessing.Value('i', 0)

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

RD = os.path.abspath("../bin/rd") + " --msa {msa} --tree {tree} --states 4 --seed {seed} --fast --silent --force"
model_file = "subst.model"
freqs_file = "freqs.model"

RUN_TEMPLATE = "run_{run_iter:0{leading_zeroes}}"
TOTAL_ITERS = 4

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
        self._params = numpy.random.rand(4,4) + 1e-2
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
        self._params += .001
        self._params /= numpy.linalg.norm(self._params, 1)
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
    def __init__(self, root_path, run_iter, trees, aligns, seed=None):
        self._run_iter = run_iter
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

        self._tree_names = []
        tree_name_counter = 0
        for tree in trees:
            if type(tree) == int:
                t = ete3.Tree()
                t.populate(tree)
                for n in t.traverse():
                    n.dist = numpy.random.exponential(0.1) + 0.005
                with open(os.path.join(self._run_path, str(tree)+".tree"), 'w') as tree_file:
                    tree_file.write(t.write(format=5))
                self._tree_names.append(str(tree))
            elif type(tree) == ete3.Tree:
                tree_name = base26_encode(tree_name_counter, len(trees))
                tree_name_counter += 1
                tree_filename = os.path.join(self._run_path, str(tree_name) + ".tree")
                with open(tree_filename,'w') as tree_file:
                    tree_file.write(tree.write(format=5))
                self._tree_names.append(tree_name)

        self._site_steps = []
        self._aligns = []
        align_name_counter = 0
        for align in aligns:
            if type(align) == int:
                self._site_steps.append(align)
            elif type(align) == list:
                align_name = base26_encode(align_name_counter, len(aligns))
                align_name_counter += 1
                self._aligns.append((align_name, align))

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

    def gen_indel_alignment(self, sites, freqs, subst, tree_filename):
        with open(tree_filename) as treef:
            tree = treef.read()
        with open('control.txt', 'w') as control_txt:
            control_txt.write(self.make_indel_control_file(freqs, subst, tree,
                sites))
        if not self.check_done_indel('.'):
            subprocess.run("indelible", stdout=subprocess.DEVNULL)
            self.set_indel_done('.')

    def run_exp(self, tree_filename, msa):
        if not self.check_done_rd('.'):
            rd_output = subprocess.run(RD.format(msa=msa,
                tree=os.path.join("../", tree_filename),
                seed=self._seed).split(' '),
                stdout=subprocess.PIPE)
            with open('rd_output_all', 'w') as logfile:
                logfile.write(rd_output.stdout.decode('utf-8'))
            lh, tree, _= rd_output.stdout.decode('utf-8').split('\n')
            with open('rd_output', 'w') as rd_outfile:
                rd_outfile.write(tree)
            with open('rd_output_lh', 'w') as rd_outfile:
                rd_outfile.write(lh)
            self.set_rd_done('.')

    def run_all(self):
        old_dir = os.getcwd()
        os.chdir(self._run_path)

        freqs, subst = self.get_model_params()
        with open('subst.model', 'w') as model_file:
            model_file.write(subst.rd_repr())

        with open('freqs.model', 'w') as model_file:
            model_file.write(freqs.rd_repr())

        for tree_name in self._tree_names:
            all_trees = []
            all_lh = [0.0]
            tree_file = os.path.join(self._run_path, str(tree_name) + ".tree")
            with open(tree_file) as tf:
                all_trees.append(tf.read()+'\n')
            for sites in self._site_steps:
                exp_dir = "{taxa}tree_{sites}sites".format(taxa=tree_name,
                            sites=sites)
                if not os.path.exists(exp_dir):
                    os.mkdir(exp_dir)
                with directory_guard(exp_dir):
                    self.gen_indel_alignment(sites, freqs, subst, tree_file)
                    self.run_exp(tree_file, 'seqs_TRUE.phy')
                    with open('rd_output') as tf:
                        all_trees.append(tf.readline())
                    with open('rd_output_lh') as tf:
                        all_lh.append(tf.readline())
            for align_name, align in self._aligns:
                exp_dir = "{taxa}tree_{align_name}align".format(taxa=tree_name,
                        align_name=align_name)
                if not os.path.exists(exp_dir):
                    os.mkdir(exp_dir)
                with directory_guard(exp_dir):
                    align_filename = str(align_name) + ".fasta"
                    with open(align_filename,'w') as align_file:
                        SeqIO.write(align, align_file, 'fasta')
                    self.run_exp(tree_file, msa=align_filename)
                    with open('rd_output') as tf:
                        all_trees.append(tf.readline())
                    with open('rd_output_lh') as tf:
                        all_lh.append(tf.readline())
            result_tree_file = "result_trees_{}_tree".format(tree_name)
            with open(result_tree_file, 'w') as rt:
               rt.write(''.join(all_trees))
            result_tree_lh = "result_trees_{}_lh".format(tree_name)
            with open(result_tree_lh, 'w') as rt:
               rt.write('\n'.join([str(f) for f in all_lh]))

            try:
                parsed_trees = [ete3.Tree(t) for t in all_trees]
            except:
                print(all_trees)
                sys.exit(1)
            true_tree = parsed_trees[0]
            self._result_trees = parsed_trees[1:]

            rfdists = []
            for i in range(1, len(parsed_trees)):
                rfdists.append(compute_root_distance(true_tree,
                    parsed_trees[i]))
            with open('rfdists_{taxa}_tree'.format(taxa=tree_name), 'w')\
                    as rf_outfile:
                rf_outfile.write(','.join([str(i) for i in self._site_steps] +
                    [a for a, _ in self._aligns]))
                rf_outfile.write('\n')
                rf_outfile.write(','.join([str(f) for f in rfdists]))
                rf_outfile.write('\n')
        PROGRESS_BAR.update(PROGRESS_BAR_ITER.value)
        PROGRESS_BAR_ITER.value += 1
        os.chdir(old_dir)

    def tree_names(self):
        return self._tree_names

    def result_trees(self):
        return self._result_trees

    def align_names(self):
        return [a for a, _ in self._aligns]

    def site_steps(self):
        return [str(s) for s in self._site_steps]

def compute_root_distance(t1, t2):
    result = t1.robinson_foulds(t2)
    return result[0]/2

def compute_average_distance(tree_names, aligns):
    for tree in tree_names:
        leading_zeroes = math.ceil(math.log10(TOTAL_ITERS))
        nrows = len(SITE_STEPS)
        totals = numpy.zeros(nrows)
        for i in range(TOTAL_ITERS):
            result_tree_file = os.path.join(RUN_TEMPLATE.format(run_iter=i,
                leading_zeroes = leading_zeroes),
                    "rfdists_{taxa}_tree".format(taxa=tree))
            with open(result_tree_file) as rfdist_file:
                rfdists = csv.DictReader(rfdist_file)
                for row in rfdists:
                    for i in range(len(SITE_STEPS)):
                        totals[i] += float(row[str(SITE_STEPS[i])])
        totals /= TOTAL_ITERS
        with open('summary_{}_tree'.format(tree), 'w') as outfile:
            outfile.write(','.join([str(i) for i in SITE_STEPS]))
            outfile.write('\n')
            outfile.write(','.join([str(f) for f in totals]))
            outfile.write('\n')

def tree_map(tree_names, trees):
    with open('tree_map', 'w') as outfile:
        for tn, t in zip(tree_names, trees):
            if type(t) != ete3.Tree:
                continue
            outfile.write(tn + ": " + t.write() + "\n")

def get_root_clade(tree):
    return sorted([n.name for n in tree.get_tree_root().children[0].traverse() if 
            n.name != ''])

def map_root_onto_main(tree_names, trees, site_steps, aligns):
    leading_zeroes = math.ceil(math.log10(TOTAL_ITERS))
    for tn, true_tree in zip(tree_names, trees):
        if type(true_tree) != ete3.Tree:
            continue
        true_tree.get_tree_root().add_features(true_root=True)
        for sites in site_steps:
            for n in true_tree.traverse():
                n.add_features(root_placement = 0)
            for i in range(TOTAL_ITERS):
                result_tree_file = os.path.join(RUN_TEMPLATE.format(leading_zeroes =
                    leading_zeroes, run_iter=i),
                 "{taxa}tree_{sites}sites".format(taxa=tn, sites=sites),
                 "rd_output")
                with open(result_tree_file) as infile:
                    result_tree = ete3.Tree(infile.readline())
                clade = get_root_clade(result_tree)
                if len(clade) == 1:
                    (true_tree & clade[0]).root_placement+=1
                else:
                    true_tree.get_common_ancestor(clade).root_placement+= 1;
            with open("{tree_name}tree_{sites}sites_mapped_tree".format(
                tree_name = tn, sites=sites), 'w') as outfile:
                outfile.write(true_tree.write(format=9,
                    features=['root_placement', 'true_root']))
        for align in aligns:
            for n in true_tree.traverse():
                n.add_features(root_placement = 0)
            for i in range(TOTAL_ITERS):
                result_tree_file = os.path.join(RUN_TEMPLATE.format(leading_zeroes =
                    leading_zeroes, run_iter=i),
                 "{taxa}tree_{align}align".format(taxa=tn, align=align),
                 "rd_output")
                with open(result_tree_file) as infile:
                    result_tree = ete3.Tree(infile.readline())
                clade = get_root_clade(result_tree)
                if len(clade) == 1:
                    (true_tree & clade[0]).root_placement+=1
                else:
                    true_tree.get_common_ancestor(clade).root_placement+= 1;
            with open("{tree_name}tree_{align}align_mapped_tree".format(
                tree_name = tn, align=align), 'w') as outfile:
                outfile.write(true_tree.write(format=9,
                    features=['root_placement', 'true_root']))

def summarize_results(path, tree_names, trees, site_steps, aligns):
    with directory_guard(path):
        tree_map(tree_names, trees)
        #compute_average_distance(tree_names)
        map_root_onto_main(tree_names, trees, site_steps, aligns)


def base26_encode(index, maximum):
    if index == 0:
        return 'a'
    iters = math.ceil(math.log(maximum, 26))
    bases = [string.ascii_lowercase[(index % ( 26 ** (e +1))) // (26 ** e)] for
        e in range(iters)]
    bases.reverse()
    return ''.join(bases)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, help='Path to store the exp',
            required=True)
    parser.add_argument('--msa', nargs='+', type=str, required=True)
    parser.add_argument('--trees', nargs='+', type=str, required=True)
    parser.add_argument('--iters', type=int, required=True)
    parser.add_argument('--procs', type=int, default=None)
    args = parser.parse_args()

    trees = []
    for tree in args.trees:
        try:
            trees.append(int(tree))
        except ValueError:
            with open(tree) as tree_file:
                trees.extend([ete3.Tree(s) for s in tree_file])

    aligns = []
    for align in args.msa:
        try:
            if not shutil.which("indelible"):
                print("Please add indelible to your path")
                sys.exit()
            aligns.append(int(align))
        except ValueError:
            aligns.append(list(SeqIO.parse(align, os.path.splitext(align)[1].strip('.'))))

    exp_path = os.path.abspath(args.path)
    TOTAL_ITERS = args.iters

    PROGRESS_BAR = progressbar.ProgressBar(max_value = TOTAL_ITERS)

    if not os.path.exists(exp_path):
        os.mkdir(exp_path)

    with directory_guard(exp_path):
        experiments = []
        for i in range(TOTAL_ITERS):
            experiments.append(exp('.', i, trees, aligns))

        PROGRESS_BAR.update(PROGRESS_BAR_ITER.value)
        PROGRESS_BAR_ITER.value += 1
        with multiprocessing.Pool(args.procs) as tp:
            tp.map(exp.run_all, experiments)
        summarize_results('.', experiments[0].tree_names(), trees,
                experiments[0].site_steps(), experiments[0].align_names())
