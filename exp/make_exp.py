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

CONTROL_FILE = """
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

RD = os.path.abspath(
    "../bin/rd"
) + " --msa {msa} --tree {tree} --states 4 --seed {seed} --silent --force"
IQTREE = "iqtree -m 12.12 -s {msa} -g {tree}"
model_file = "subst.model"
freqs_file = "freqs.model"

RUN_TEMPLATE = "run_{run_iter:0{leading_zeroes}}"
TOTAL_ITERS = 0


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
        self._params = numpy.random.rand(4, 4) + 1e-2
        self._params -= numpy.diag(numpy.diag(self._params))
        self._params -= numpy.diagflat(
            numpy.dot(self._params, numpy.ones((4, 1))))

    def indel_repr(self):
        p = []
        for i in range(4):
            for j in range(4):
                if i == j:
                    continue
                p.append(self._params[i][j])
        return ' '.join([str(f) for f in p])

    def rd_repr(self):
        P = numpy.array([[0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0],
                         [0, 0, 1, 0]])
        tmp = numpy.inner(numpy.inner(P.T, self._params), P)
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
        P = numpy.array([[0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0],
                         [0, 0, 1, 0]])
        p = numpy.inner(self._params, P)
        return ','.join([str(f) for f in p])


class exp:
    def __init__(self, root_path, run_iter, trees, aligns, seed=None):
        self._run_iter = run_iter
        leading_zeroes = math.ceil(math.log10(TOTAL_ITERS))
        self._run_path = os.path.abspath(
            os.path.join(
                root_path,
                RUN_TEMPLATE.format(run_iter=run_iter,
                                    leading_zeroes=leading_zeroes)))
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
                t.unroot()
                for n in t.traverse():
                    n.dist = numpy.random.exponential(0.1) + 0.005
                with open(os.path.join(self._run_path,
                                       str(tree) + ".tree"), 'w') as tree_file:
                    tree_file.write(t.write(format=5))
                self._tree_names.append(str(tree))
            elif type(tree) == ete3.Tree:
                tree_name = base26_encode(tree_name_counter, len(trees))
                tree_name_counter += 1
                tree_filename = os.path.join(self._run_path,
                                             str(tree_name) + ".tree")
                unrooted_tree = tree.copy()
                unrooted_tree.unroot()
                with open(tree_filename, 'w') as tree_file:
                    tree_file.write(unrooted_tree.write(format=5))
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
    def set_iqtree_done(path):
        with open(os.path.join(path, ".done"), 'a') as done_file:
            done_file.write("iqtree:" + datetime.datetime.now().isoformat())

    @staticmethod
    def check_done_rd(path):
        if os.path.exists(os.path.join(path, '.done')):
            with open(os.path.join(path, ".done")) as done_file:
                for line in done_file:
                    if line.find('rd') != -1:
                        return True
        return False

    @staticmethod
    def check_done_iqtree(path):
        if os.path.exists(os.path.join(path, '.done')):
            with open(os.path.join(path, ".done")) as done_file:
                for line in done_file:
                    if line.find('iqtree') != -1:
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
        return CONTROL_FILE.format(freq_params=freqs.indel_repr(),
                                   model_params=subst.indel_repr(),
                                   tree=tree,
                                   sites=sites,
                                   randomseed=self._seed)

    def gen_indel_alignment(self, sites, freqs, subst, tree_filename):
        with open(tree_filename) as treef:
            tree = treef.read()
        with open('control.txt', 'w') as control_txt:
            control_txt.write(
                self.make_indel_control_file(freqs, subst, tree, sites))
        if not self.check_done_indel('.'):
            subprocess.run("indelible", stdout=subprocess.DEVNULL)
            self.set_indel_done('.')

    def run_rd(self, tree_filename, msa):
        rd_output = subprocess.run(RD.format(msa=msa,
                                             tree=os.path.join(
                                                 "../", tree_filename),
                                             seed=self._seed).split(' '),
                                   stdout=subprocess.PIPE)
        with open('rd_output_all', 'w') as logfile:
            logfile.write(rd_output.stdout.decode('utf-8'))
        lh, tree, _ = rd_output.stdout.decode('utf-8').split('\n')
        with open('rd_output', 'w') as rd_outfile:
            rd_outfile.write(tree)
        with open('rd_output_lh', 'w') as rd_outfile:
            rd_outfile.write(lh)
        self.set_rd_done('.')

    def run_iqtree(self, tree_filename, msa):
        subprocess.run(IQTREE.format(msa=msa, tree=tree_filename).split(),
                       stdout=subprocess.DEVNULL)
        self.set_iqtree_done('.')

    def run_exp(self, tree_filename, msa):
        if not self.check_done_rd('.'):
            self.run_rd(tree_filename, msa)
        if not self.check_done_iqtree('.'):
            self.run_iqtree(tree_filename, msa)

    def run_all(self):
        old_dir = os.getcwd()
        os.chdir(self._run_path)

        freqs, subst = self.get_model_params()
        with open('subst.model', 'w') as model_file:
            model_file.write(subst.rd_repr())

        with open('freqs.model', 'w') as model_file:
            model_file.write(freqs.rd_repr())

        for tree_name in self._tree_names:
            all_trees_rd = []
            all_trees_iqtree = []
            all_lh_rd = [0.0]
            tree_file = os.path.join(self._run_path, str(tree_name) + ".tree")
            with open(tree_file) as tf:
                true_tree_newick = tf.read()
            for sites in self._site_steps:
                exp_dir = "{taxa}tree_{sites}sites".format(taxa=tree_name,
                                                           sites=sites)
                if not os.path.exists(exp_dir):
                    os.mkdir(exp_dir)
                with directory_guard(exp_dir):
                    self.gen_indel_alignment(sites, freqs, subst, tree_file)
                    self.run_exp(tree_file, 'seqs_TRUE.phy')
                    with open('rd_output') as tf:
                        all_trees_rd.append(tf.readline())
                    with open('rd_output_lh') as tf:
                        all_lh_rd.append(tf.readline())
            for align_name, align in self._aligns:
                exp_dir = "{taxa}tree_{align_name}align".format(
                    taxa=tree_name, align_name=align_name)
                if not os.path.exists(exp_dir):
                    os.mkdir(exp_dir)
                with directory_guard(exp_dir):
                    align_filename = str(align_name) + ".fasta"
                    with open(align_filename, 'w') as align_file:
                        SeqIO.write(align, align_file, 'fasta')
                    self.run_exp(tree_file, msa=align_filename)
                    with open('rd_output') as tf:
                        all_trees_rd.append(tf.readline())
                    with open('rd_output_lh') as tf:
                        all_lh_rd.append(tf.readline())
                    with open('{}.treefile'.format(
                            align_filename)) as iqtree_treefile:
                        all_trees_iqtree.append(iqtree_treefile.read())
            result_tree_file_rd = "result_trees_{}_tree_rd".format(tree_name)
            result_tree_file_iqtree = "result_trees_{}_tree_iqtree".format(
                tree_name)
            with open(result_tree_file_rd, 'w') as rt:
                rt.write(''.join(all_trees_rd))
            with open(result_tree_file_iqtree, 'w') as rt:
                rt.write(''.join(all_trees_iqtree))
            result_tree_lh = "result_trees_{}_lh_rd".format(tree_name)
            with open(result_tree_lh, 'w') as rt:
                rt.write('\n'.join([str(f) for f in all_lh_rd]))

            try:
                parsed_trees = [ete3.Tree(t) for t in all_trees_rd]
            except:
                print(all_trees_rd)
                sys.exit(1)
            true_tree = ete3.Tree(true_tree_newick)
            self._result_trees = parsed_trees

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


def tree_map(tree_names, trees):
    with open('tree_map', 'w') as outfile:
        for tn, t in zip(tree_names, trees):
            if type(t) != ete3.Tree:
                continue
            outfile.write(tn + ": " + t.write() + "\n")


def get_left_clade(tree):
    return sorted([
        n.name for n in tree.get_tree_root().children[0].traverse()
        if n.name != ''
    ])


def get_right_clade(tree):
    return sorted([
        n.name for n in tree.get_tree_root().children[1].traverse()
        if n.name != ''
    ])


def get_root_clades(tree):
    left_clade = get_left_clade(tree)
    right_clade = get_right_clade(tree)
    return (left_clade, right_clade)


def extract_node_with_clade(tree, clade):
    if len(clade) == 1:
        return (tree & clade[0])
    return tree.get_common_ancestor(clade)


def get_mapped_node(true_tree, inferred_tree):
    left_clade, right_clade = get_root_clades(inferred_tree)

    left_node = extract_node_with_clade(true_tree, left_clade)
    right_node = extract_node_with_clade(true_tree, right_clade)

    if left_node == right_node:
        return left_node

    if left_node in right_node.get_descendants():
        return left_node

    if right_node in left_node.get_descendants():
        return right_node

    if left_node.is_leaf():
        return left_node
    if right_node.is_leaf():
        return right_node

    return true_tree.get_tree_root()


def get_root_distance_metric(true_tree, inferred_tree):
    common_node_tt, _ = true_tree.get_closest_leaf()
    common_node_it = inferred_tree & common_node_tt.name
    cn_tt_dist = true_tree.get_distance(common_node_tt,
                                        true_tree.get_tree_root())
    cn_it_dist = inferred_tree.get_distance(common_node_it,
                                            inferred_tree.get_tree_root())
    return numpy.abs(cn_tt_dist - cn_it_dist)


def get_root_distance_toplogical(true_tree, inferred_tree):
    common_node_tt, _ = true_tree.get_closest_leaf()
    common_node_it = inferred_tree & common_node_tt.name
    cn_tt_dist = true_tree.get_distance(common_node_tt,
                                        true_tree.get_tree_root(),
                                        topology_only=True)
    cn_it_dist = inferred_tree.get_distance(common_node_it,
                                            inferred_tree.get_tree_root(),
                                            topology_only=True)
    return numpy.abs(cn_tt_dist - cn_it_dist)


def compute_distances(tree_names, trees, site_steps, aligns):
    leading_zeroes = math.ceil(math.log10(TOTAL_ITERS))
    for tn, true_tree in zip(tree_names, trees):
        for sites in site_steps:
            rf_topo_iqtree = []
            rf_topo_rd = []
            rf_metric_iqtree = []
            rf_metric_rd = []
            for i in range(TOTAL_ITERS):
                if type(true_tree) != ete3.Tree:
                    with open(
                            os.path.join(
                                RUN_TEMPLATE.format(
                                    leading_zeroes=leading_zeroes, run_iter=i),
                                "{}.tree".format(tn))) as infile:
                        true_tree = ete3.Tree(infile.read())
                result_tree_file_rd = os.path.join(
                    RUN_TEMPLATE.format(leading_zeroes=leading_zeroes,
                                        run_iter=i),
                    "{taxa}tree_{sites}sites".format(taxa=tn,
                                                     sites=sites), "rd_output")
                result_tree_file_iqtree = os.path.join(
                    RUN_TEMPLATE.format(leading_zeroes=leading_zeroes,
                                        run_iter=i),
                    "{taxa}tree_{sites}sites".format(taxa=tn, sites=sites),
                    "seqs_TRUE.phy.treefile")
                with open(result_tree_file_rd) as infile:
                    result_tree_rd = ete3.Tree(infile.readline())
                with open(result_tree_file_iqtree) as infile:
                    result_tree_iqtree = ete3.Tree(infile.readline())

                rf_topo_rd.append(
                    get_root_distance_toplogical(true_tree, result_tree_rd))
                rf_metric_rd.append(
                    get_root_distance_metric(true_tree, result_tree_rd))
                rf_topo_iqtree.append(
                    get_root_distance_toplogical(true_tree,
                                                 result_tree_iqtree))
                rf_metric_iqtree.append(
                    get_root_distance_metric(true_tree, result_tree_iqtree))

            tree_size = len(true_tree.get_leaves())
            with open(
                    "{tree_name}tree_{sites}sites_rf_dists".format(
                        tree_name=tn, sites=sites), 'w') as outfile:
                outfile.write('method,topo,ntopo,metric\n')
                outfile.write('rd,{},{},{}\n'.format(
                    numpy.mean(rf_topo_rd),
                    numpy.mean(rf_topo_rd) / tree_size,
                    numpy.mean(rf_metric_rd)))
                outfile.write('iqtree,{},{},{}\n'.format(
                    numpy.mean(rf_topo_iqtree),
                    numpy.mean(rf_topo_iqtree) / tree_size,
                    numpy.mean(rf_metric_rd)))
            with open(
                    "{tree_name}tree_{sites}sites_rf_hist".format(tree_name=tn,
                                                                  sites=sites),
                    'w') as outfile:
                outfile.write('iqtree_n_dist,rd_n_dist\n')
                for rd, iq in zip(rf_topo_rd, rf_topo_iqtree):
                    outfile.write("{},{},{},{}\n".format(
                        iq, iq / tree_size, rd, rd / tree_size))

        if type(true_tree) != ete3.Tree:
            with open(
                    os.path.join(
                        RUN_TEMPLATE.format(leading_zeroes=leading_zeroes,
                                            run_iter=i),
                        "{}.tree".format(tn))) as infile:
                true_tree = ete3.Tree(infile.read())
        for align in aligns:
            rf_topo_iqtree = []
            rf_topo_rd = []
            rf_metric_iqtree = []
            rf_metric_rd = []
            for i in range(TOTAL_ITERS):
                result_tree_file_rd = os.path.join(
                    RUN_TEMPLATE.format(leading_zeroes=leading_zeroes,
                                        run_iter=i),
                    "{taxa}tree_{align}align".format(taxa=tn,
                                                     align=align), "rd_output")
                result_tree_file_iqtree = os.path.join(
                    RUN_TEMPLATE.format(leading_zeroes=leading_zeroes,
                                        run_iter=i),
                    "{taxa}tree_{align}align".format(taxa=tn, align=align),
                    "{align}.fasta.treefile".format(align=align))
                with open(result_tree_file_rd) as infile:
                    result_tree_rd = ete3.Tree(infile.readline())
                with open(result_tree_file_iqtree) as infile:
                    result_tree_iqtree = ete3.Tree(infile.readline())
                rf_topo_rd.append(
                    get_root_distance_toplogical(true_tree, result_tree_rd))
                rf_metric_rd.append(
                    get_root_distance_metric(true_tree, result_tree_rd))
                rf_topo_iqtree.append(
                    get_root_distance_toplogical(true_tree,
                                                 result_tree_iqtree))
                rf_metric_iqtree.append(
                    get_root_distance_metric(true_tree, result_tree_iqtree))

            tree_size = len(true_tree.get_leaves())
            with open(
                    "{tree_name}tree_{align}align_rf_dists".format(
                        tree_name=tn, align=align), 'w') as outfile:
                outfile.write('method,topo,ntopo,metric\n')
                outfile.write('rd,{},{},{}\n'.format(
                    numpy.mean(rf_topo_rd),
                    numpy.mean(rf_topo_rd) / tree_size,
                    numpy.mean(rf_metric_rd)))
                outfile.write('iqtree,{},{},{}\n'.format(
                    numpy.mean(rf_topo_iqtree),
                    numpy.mean(rf_topo_iqtree) / tree_size,
                    numpy.mean(rf_metric_rd)))
            with open(
                    "{tree_name}tree_{align}align_rf_hist".format(tree_name=tn,
                                                                  align=align),
                    'w') as outfile:
                outfile.write('iqtree_dist,iqtree_n_dist,rd_dist,rd_n_dist\n')
                for rd, iq in zip(rf_topo_rd, rf_topo_iqtree):
                    outfile.write("{},{},{},{}\n".format(
                        iq, iq / tree_size, rd, rd / tree_size))


def map_root_onto_main(tree_names, trees, site_steps, aligns):
    leading_zeroes = math.ceil(math.log10(TOTAL_ITERS))
    for tn, true_tree in zip(tree_names, trees):
        if type(true_tree) != ete3.Tree:
            continue
        for sites in site_steps:
            for n in true_tree.traverse():
                n.add_features(root_placement_rd=0)
                n.add_features(root_placement_iqtree=0)
            for i in range(TOTAL_ITERS):
                result_tree_file_rd = os.path.join(
                    RUN_TEMPLATE.format(leading_zeroes=leading_zeroes,
                                        run_iter=i),
                    "{taxa}tree_{sites}sites".format(taxa=tn,
                                                     sites=sites), "rd_output")
                result_tree_file_iqtree = os.path.join(
                    RUN_TEMPLATE.format(leading_zeroes=leading_zeroes,
                                        run_iter=i),
                    "{taxa}tree_{sites}sites".format(taxa=tn, sites=sites),
                    "seqs_TRUE.phy.treefile")
                with open(result_tree_file_rd) as infile:
                    result_tree_rd = ete3.Tree(infile.readline())
                with open(result_tree_file_iqtree) as infile:
                    result_tree_iqtree = ete3.Tree(infile.readline())
                clade_rd = get_mapped_node(true_tree, result_tree_rd)
                clade_iqtree = get_mapped_node(true_tree, result_tree_iqtree)

                clade_rd.root_placement_rd += 1
                clade_iqtree.root_placement_iqtree += 1
            with open(
                    "{tree_name}tree_{sites}sites_mapped_tree".format(
                        tree_name=tn, sites=sites), 'w') as outfile:
                outfile.write(
                    true_tree.write(format=9,
                                    features=[
                                        'root_placement_rd',
                                        'root_placement_iqtree'
                                    ]))

        for align in aligns:
            for n in true_tree.traverse():
                n.add_features(root_placement_rd=0)
                n.add_features(root_placement_iqtree=0)
            for i in range(TOTAL_ITERS):
                result_tree_file_rd = os.path.join(
                    RUN_TEMPLATE.format(leading_zeroes=leading_zeroes,
                                        run_iter=i),
                    "{taxa}tree_{align}align".format(taxa=tn,
                                                     align=align), "rd_output")
                result_tree_file_iqtree = os.path.join(
                    RUN_TEMPLATE.format(leading_zeroes=leading_zeroes,
                                        run_iter=i),
                    "{taxa}tree_{align}align".format(taxa=tn, align=align),
                    "{align}.fasta.treefile".format(align=align))
                with open(result_tree_file_rd) as infile:
                    result_tree_rd = ete3.Tree(infile.readline())
                with open(result_tree_file_iqtree) as infile:
                    result_tree_iqtree = ete3.Tree(infile.readline())
                clade_rd = get_mapped_node(true_tree, result_tree_rd)
                clade_iqtree = get_mapped_node(true_tree, result_tree_iqtree)

                clade_rd.root_placement_rd += 1
                clade_iqtree.root_placement_iqtree += 1
            with open(
                    "{tree_name}tree_{align}align_mapped_tree".format(
                        tree_name=tn, align=align), 'w') as outfile:
                outfile.write(
                    true_tree.write(format=9,
                                    features=[
                                        'root_placement_rd',
                                        'root_placement_iqtree'
                                    ]))


def summarize_results(path, tree_names, trees, site_steps, aligns):
    with directory_guard(path):
        tree_map(tree_names, trees)
        map_root_onto_main(tree_names, trees, site_steps, aligns)
        compute_distances(tree_names, trees, site_steps, aligns)


def base26_encode(index, maximum):
    if index == 0:
        return 'a'
    iters = math.ceil(math.log(maximum, 26))
    bases = [
        string.ascii_lowercase[(index % (26**(e + 1))) // (26**e)]
        for e in range(iters)
    ]
    bases.reverse()
    return ''.join(bases)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--path',
                        type=str,
                        help='Path to store the exp',
                        required=True)
    parser.add_argument('--msa', nargs='+', type=str, required=True)
    parser.add_argument('--trees', nargs='+', type=str, required=True)
    parser.add_argument('--iters', type=int, required=True)
    parser.add_argument('--procs', type=int, default=None)
    args = parser.parse_args()

    if not shutil.which("iqtree"):
        print("Please add iqtree to your path")
        sys.exit()

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
            aligns.append(
                list(SeqIO.parse(align,
                                 os.path.splitext(align)[1].strip('.'))))

    exp_path = os.path.abspath(args.path)
    TOTAL_ITERS = args.iters

    PROGRESS_BAR = progressbar.ProgressBar(max_value=TOTAL_ITERS)

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
                          experiments[0].site_steps(),
                          experiments[0].align_names())
