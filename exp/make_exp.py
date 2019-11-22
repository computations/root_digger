#!/usr/bin/env python3

import argparse
import csv
import datetime
import math
import multiprocessing
from multiprocessing.pool import ThreadPool
import os
import random
import shutil
import string
import subprocess
import sys
import json
import itertools

import Bio
import ete3
import numpy
import progressbar
from Bio import SeqIO

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
) + " --msa {msa} --tree {tree} --states 4 --seed {seed} --force"
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
    def __init__(self,
                 root_path,
                 run_iter,
                 trees,
                 aligns,
                 run_rd=True,
                 run_iq=True,
                 seed=None):
        self._run_rd = run_rd
        self._run_iq = run_iq
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
                for n in t.traverse():
                    n.dist = numpy.random.exponential(0.1) + 0.005
                with open(os.path.join(self._run_path,
                                       str(tree) + ".rtree"),
                          'w') as tree_file:
                    tree_file.write(t.write(format=5))
                t.unroot()
                with open(os.path.join(self._run_path,
                                       str(tree) + ".tree"), 'w') as tree_file:
                    tree_file.write(t.write(format=5))
                self._tree_names.append(str(tree))
            elif type(tree) == ete3.Tree:
                tree_name = base26_encode(tree_name_counter, len(trees))
                tree_name_counter += 1
                tree_filename = os.path.join(self._run_path,
                                             str(tree_name) + ".tree")
                rtree_filename = os.path.join(self._run_path,
                                              str(tree_name) + ".rtree")
                unrooted_tree = tree.copy()
                with open(rtree_filename, 'w') as tree_file:
                    tree_file.write(unrooted_tree.write(format=5))
                unrooted_tree.unroot()
                with open(tree_filename, 'w') as tree_file:
                    tree_file.write(unrooted_tree.write(format=5))
                self._tree_names.append(tree_name)

        self._site_steps = []
        self._aligns = []
        self._exp_keys = set()
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
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        with open('rd_output', 'w') as logfile:
            logfile.write(rd_output.stdout.decode('utf-8'))
        with open('rd_output_err', 'w') as logfile:
            logfile.write(rd_output.stderr.decode('utf-8'))
        self.set_rd_done('.')

    def run_iqtree(self, tree_filename, msa):
        subprocess.run(IQTREE.format(msa=msa, tree=tree_filename).split(),
                       stdout=subprocess.DEVNULL)
        self.set_iqtree_done('.')

    def run_exp(self, tree_filename, msa):
        if not self.check_done_rd('.') and self._run_rd:
            self.run_rd(tree_filename, msa)
        if not self.check_done_iqtree('.') and self._run_iq:
            self.run_iqtree(tree_filename, msa)

    def run_all(self):
        old_dir = os.getcwd()
        os.chdir(self._run_path)
        self._rd_results = {}
        self._iqtree_results = {}

        freqs, subst = self.get_model_params()
        with open('subst.model', 'w') as model_file:
            model_file.write(subst.indel_repr())

        with open('freqs.model', 'w') as model_file:
            model_file.write(freqs.indel_repr())

        for tree_name in self._tree_names:
            tree_file = os.path.join(self._run_path, str(tree_name) + ".tree")
            rtree_file = os.path.join(self._run_path,
                                      str(tree_name) + ".rtree")
            with open(rtree_file) as tf:
                true_tree_newick = tf.read()
                true_tree_ete = ete3.Tree(true_tree_newick)
            for sites in self._site_steps:
                exp_dir = "{taxa}tree_{sites}sites".format(taxa=tree_name,
                                                           sites=sites)
                if not os.path.exists(exp_dir):
                    os.mkdir(exp_dir)
                with directory_guard(exp_dir):
                    self.gen_indel_alignment(sites, freqs, subst, rtree_file)
                    self.run_exp(tree_file, 'seqs_TRUE.phy')
                exp_key = (tree_name, sites)
                if self._run_rd:
                    self._rd_results[exp_key] = rd_result(
                        exp_dir, true_tree_ete)
                if self._run_iq:
                    self._iqtree_results[exp_key] = iqtree_result(
                        exp_dir, true_tree_ete, 'seqs_TRUE.phy')
                if exp_key not in self._exp_keys:
                    self._exp_keys.add(exp_key)

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
                exp_key = (tree_name, align_name)
                if self._run_rd:
                    self._rd_results[exp_key] = rd_result(
                        exp_dir, true_tree_ete)
                if self._run_iq:
                    self._iqtree_results[exp_key] = iqtree_result(
                        exp_dir, true_tree_ete, align_filename)
                if exp_key not in self._exp_keys:
                    self._exp_keys.add(exp_key)

        PROGRESS_BAR.update(PROGRESS_BAR_ITER.value)
        PROGRESS_BAR_ITER.value += 1
        os.chdir(old_dir)
        return self

    def tree_names(self):
        return self._tree_names

    def result_trees(self):
        return self._result_trees

    def align_names(self):
        return [a for a, _ in self._aligns]

    def site_steps(self):
        return [str(s) for s in self._site_steps]

    def rd_results(self):
        if self._run_rd:
            return self._rd_results
        return None

    def iqtree_results(self):
        if self._run_iq:
            return self._iqtree_results
        return None

    def exp_keys(self):
        return self._exp_keys


class result:
    @property
    def time(self):
        return self._time

    @property
    def tree(self):
        return self._tree.write(format=5)

    @property
    def lh(self):
        return self._final_lh

    @property
    def root_distance(self):
        return self._root_distance

    @property
    def normalized_root_distance(self):
        return self._normalized_root_distance

    @property
    def path_distance(self):
        return self._path_distance

    @property
    def normalized_path_distance(self):
        return self._normalized_path_distance

    @property
    def right_clade_list(self):
        return self.get_right_clade_list(self._tree)

    @property
    def left_clade_list(self):
        return self.get_left_clade_list(self._tree)

    @property
    def right_clade(self):
        return self._right_clade

    @property
    def left_clade(self):
        return self._left_clade

    @property
    def true_tree(self):
        return self._true_tree

    def _calculate_distances(self, true_tree, inferred_tree):
        tree_size = len(true_tree.get_leaves())
        self._root_distance = get_root_distance_toplogical(
            true_tree, self._tree)
        self._normalized_root_distance = self._root_distance / tree_size
        self._path_distance = get_root_distance_metric(true_tree, self._tree)
        self._normalized_path_distance = self._path_distance / tree_size

    def get(self):
        return {
            'time': self.time,
            'tree': self.tree,
            'lh': self.lh,
            'root_distance': self.root_distance,
            'path_distance': self.path_distance,
            'normalized_root_distance': self.normalized_root_distance,
            'normalized_path_distance': self.normalized_path_distance,
            'right_clade': self.right_clade,
            'left_clade': self.left_clade,
        }

    @staticmethod
    def get_left_clade_list(tree):
        return sorted([
            n.name for n in tree.get_tree_root().children[0].traverse()
            if n.name != ''
        ])

    @staticmethod
    def get_right_clade_list(tree):
        return sorted([
            n.name for n in tree.get_tree_root().children[1].traverse()
            if n.name != ''
        ])

    @staticmethod
    def get_left_clade(tree):
        return ''.join(
            sorted([
                n.name for n in tree.get_tree_root().children[0].traverse()
                if n.name != ''
            ]))

    @staticmethod
    def get_right_clade(tree):
        return ''.join(
            sorted([
                n.name for n in tree.get_tree_root().children[1].traverse()
                if n.name != ''
            ]))


class rd_result(result):
    def __init__(self, directory, true_tree):
        with directory_guard(directory):
            with open('rd_output') as results_file:
                results_string = results_file.read().split('\n')
        self._time = rd_result._read_time(results_string)
        self._tree = rd_result._read_tree(results_string)
        self._final_lh = rd_result._read_lh(results_string)
        self._calculate_distances(true_tree, self._tree)
        self._right_clade = result.get_right_clade(self._tree)
        self._left_clade = result.get_left_clade(self._tree)
        self._true_tree = true_tree

    @staticmethod
    def _read_time(results):
        time_line = results[-2]
        start_index = len('Inference took: ')
        end_index = start_index + time_line[start_index:].find('s')
        return float(time_line[start_index:end_index])

    @staticmethod
    def _read_tree(results):
        tree_string = results[-3]
        return ete3.Tree(tree_string)

    @staticmethod
    def _read_lh(results):
        lh_string = results[-4]
        start_index = lh_string.find(':') + 1
        return float(lh_string[start_index:])


class iqtree_result(result):
    def __init__(self, directory, true_tree, prefix):
        with directory_guard(directory):
            with open('{}.treefile'.format(prefix)) as tree_file:
                self._tree = iqtree_result._read_tree(tree_file.read())
            with open('{}.iqtree'.format(prefix)) as iqtree_file:
                for line in iqtree_file:
                    if 'Log-likelihood of the tree:' in line:
                        self._final_lh = iqtree_result._read_lh(line)
                    if 'Total wall-clock time used:' in line:
                        self._time = iqtree_result._read_time(line)
        self._calculate_distances(true_tree, self._tree)
        self._right_clade = result.get_right_clade(self._tree)
        self._left_clade = result.get_left_clade(self._tree)
        self._true_tree = true_tree

    @staticmethod
    def _read_tree(tree_string):
        return ete3.Tree(tree_string)

    @staticmethod
    def _read_lh(lh_string):
        start_index = len('Log-likelihood of the tree: ')
        end_index = lh_string.rfind('(') - 1
        return float(lh_string[start_index:end_index])

    @staticmethod
    def _read_time(time_string):
        start_index = len('Total wall-clock time used: ')
        end_index = start_index + time_string[start_index:].find(' ')
        return float(time_string[start_index:end_index])


class summary_row:
    def __init__(self, exp, stats):
        self._values = {}
        self._values['tree'] = exp[0]
        self._values['alignment'] = exp[1]
        for stat_name, stat_value in stats.items():
            rd_key = 'rd_' + stat_name
            iq_key = 'iq_' + stat_name
            if stat_name != 'true_tree':
                self._values[rd_key] = stat_value['rd']
                self._values[iq_key] = stat_value['iq']
            else:
                self._values['true_tree'] = stat_value.write(format=5)

    def make_row(self, header):
        return ','.join([str(s) for s in self.make_line(header)])

    def make_line(self, header):
        return [
            self._values[h] if type(self._values[h]) != dict
            and type(self._values[h]) != list else json.dumps(self._values[h])
            for h in header
        ]

    def dict(self):
        return self._values


class summary:
    def __init__(self, experiments):
        self._experiments = summary._extract_exp_keys(experiments)
        self._rd_results = [e.rd_results() for e in experiments]
        self._iq_results = [e.iqtree_results() for e in experiments]

    def make_header(self):
        return [
            'tree',
            'alignment',
            'rd_time_mean',
            'iq_time_mean',
            'rd_time_median',
            'iq_time_median',
            'rd_time_std',
            'iq_time_std',
            'rd_root_distance_mean',
            'iq_root_distance_mean',
            'rd_root_distance_median',
            'iq_root_distance_median',
            'rd_root_distance_std',
            'iq_root_distance_std',
            'rd_path_distance_mean',
            'iq_path_distance_mean',
            'rd_path_distance_median',
            'iq_path_distance_median',
            'rd_path_distance_std',
            'iq_path_distance_std',
            'rd_left_clade',
            'iq_left_clade',
            'rd_right_clade',
            'iq_right_clade',
        ]

    def write(self, prefix):
        csv_filename = prefix + '.csv'
        header = self.make_header()
        with open(csv_filename, 'w') as results_file:
            results_file.write(','.join(header))
            results_file.write('\n')
            for row in self.generate_rows():
                results_file.write(row.make_row(header))
                results_file.write('\n')
        json_filename = prefix + '.json'
        with open(json_filename, 'w') as results_file:
            json.dump([r.dict() for r in self.generate_rows(json=True)],
                      results_file)

    def generate_rows(self, json=False):
        for k in self._experiments:
            stats = {}
            stats['time_mean'] = self.mean_times(k)
            stats['time_median'] = self.median_times(k)
            stats['time_std'] = self.std_times(k)
            stats['root_distance_mean'] = self.mean_root_distance(k)
            stats['root_distance_median'] = self.median_root_distance(k)
            stats['root_distance_std'] = self.std_root_distance(k)
            stats['path_distance_mean'] = self.mean_root_distance(k)
            stats['path_distance_median'] = self.median_root_distance(k)
            stats['path_distance_std'] = self.std_root_distance(k)
            if json:
                stats['true_tree'] = self.true_tree(k)
                stats['left_clade'] = self.left_clades_list(k)
                stats['right_clade'] = self.right_clades_list(k)
            else:
                stats['left_clade'] = self.left_clades(k)
                stats['right_clade'] = self.right_clades(k)
            yield summary_row(k, stats)

    def true_tree(self, k):
        if not self._rd_results[0][k].true_tree is None:
            return self._rd_results[0][k].true_tree
        return self._iq_results[k].true_tree

    def mean_times(self, k):
        return {
            'rd': summary._mean_attr(self._rd_results, k, 'time'),
            'iq': summary._mean_attr(self._iq_results, k, 'time')
        }

    def median_times(self, k):
        return {
            'rd': summary._median_attr(self._rd_results, k, 'time'),
            'iq': summary._median_attr(self._iq_results, k, 'time')
        }

    def std_times(self, k):
        return {
            'rd': summary._stddev_attr(self._rd_results, k, 'time'),
            'iq': summary._stddev_attr(self._iq_results, k, 'time')
        }

    def mean_root_distance(self, k):
        return {
            'rd': summary._mean_attr(self._rd_results, k, 'root_distance'),
            'iq': summary._mean_attr(self._iq_results, k, 'root_distance')
        }

    def median_root_distance(self, k):
        return {
            'rd': summary._median_attr(self._rd_results, k, 'root_distance'),
            'iq': summary._median_attr(self._iq_results, k, 'root_distance')
        }

    def std_root_distance(self, k):
        return {
            'rd': summary._stddev_attr(self._rd_results, k, 'root_distance'),
            'iq': summary._stddev_attr(self._iq_results, k, 'root_distance')
        }

    def mean_path_distance(self, k):
        return {
            'rd': summary._mean_attr(self._rd_results, k, 'path_distance'),
            'iq': summary._mean_attr(self._iq_results, k, 'path_distance')
        }

    def median_path_distance(self, k):
        return {
            'rd': summary._median_attr(self._rd_results, k, 'path_distance'),
            'iq': summary._median_attr(self._iq_results, k, 'path_distance')
        }

    def std_path_distance(self, k):
        return {
            'rd': summary._stddev_attr(self._rd_results, k, 'path_distance'),
            'iq': summary._stddev_attr(self._iq_results, k, 'path_distance')
        }

    def left_clades(self, k):
        return {
            'rd': summary._count_values(self._rd_results, k, 'left_clade'),
            'iq': summary._count_values(self._iq_results, k, 'left_clade'),
        }

    def right_clades(self, k):
        return {
            'rd': summary._count_values(self._rd_results, k, 'right_clade'),
            'iq': summary._count_values(self._iq_results, k, 'right_clade'),
        }

    def left_clades_list(self, k):
        return {
            'rd': summary._create_list(self._rd_results, k, 'left_clade_list'),
            'iq': summary._create_list(self._iq_results, k, 'left_clade_list'),
        }

    def right_clades_list(self, k):
        return {
            'rd': summary._create_list(self._rd_results, k,
                                       'right_clade_list'),
            'iq': summary._create_list(self._iq_results, k,
                                       'right_clade_list'),
        }

    @staticmethod
    def _create_list(results, e, key):
        if not None in results:
            items = []
            for r in results:
                k = getattr(r[e], key)
                items.append(k)
            return items
        return None

    @staticmethod
    def _count_values(results, e, key):
        if not None in results:
            counts = {}
            for r in results:
                k = getattr(r[e], key)
                if k in counts:
                    counts[k] += 1
                else:
                    counts[k] = 1
            return counts
        return None

    @staticmethod
    def _mean_attr(results, e, key):
        if not None in results:
            return numpy.mean([getattr(r[e], key) for r in results])
        return None

    @staticmethod
    def _median_attr(results, e, key):
        if not None in results:
            return numpy.median([getattr(r[e], key) for r in results])
        return None

    @staticmethod
    def _stddev_attr(results, e, key):
        if not None in results:
            return numpy.std([getattr(r[e], key) for r in results])
        return None

    @staticmethod
    def _extract_exp_keys(experiments):
        keys = set()
        for e in experiments:
            keys = keys | e.exp_keys()
        return keys


def tree_map(tree_names, trees):
    with open('tree_map', 'w') as outfile:
        for tn, t in zip(tree_names, trees):
            if type(t) != ete3.Tree:
                continue
            outfile.write(tn + ": " + t.write() + "\n")


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


def extract_node_with_clade_pair(tree, clade_pair):
    if len(clade_pair[0]) == 1:
        return (tree & clade_pair[0][0])
    if len(clade_pair[1]) == 1:
        return (tree & clade_pair[1][0])
    clade = clade_pair[0] if len(clade_pair[0]) <= len(clade_pair[1]) else\
        clade_pair[1]
    ca0 = tree.get_common_ancestor(clade_pair[0])
    ca1 = tree.get_common_ancestor(clade_pair[1])
    if ca0 in tree.children and ca1 in tree.children:
        return tree
    return tree.get_common_ancestor(clade)


def map_root_rd(true_tree, left_clades, right_clades):
    for clade_pair in zip(left_clades, right_clades):
        node = extract_node_with_clade_pair(true_tree, clade_pair)
        if node is None:
            continue
        if not hasattr(node, 'rd_map'):
            node.add_features(rd_map=1)
        else:
            node.rd_map += 1


def map_root_iq(true_tree, left_clades, right_clades):
    for clade_pair in zip(left_clades, right_clades):
        node = extract_node_with_clade_pair(true_tree, clade_pair)
        if node is None:
            continue
        if not hasattr(node, 'iq_map'):
            node.add_features(iq_map=1)
        else:
            node.iq_map += 1


def produce_mapped_root_images(json_results_filename, map_iqtree=True,
        map_rd=True):
    with open(json_results_filename) as json_file:
        results = json.load(json_file)
    for result in results:
        mapped_tree = ete3.Tree(result['true_tree'], format=5)
        if map_rd:
            map_root_rd(mapped_tree, result['rd_left_clade'],
                        result['rd_right_clade'])
        if map_iqtree:
            map_root_iq(mapped_tree, result['iq_left_clade'],
                        result['iq_right_clade'])

        def layout(node):
            if map_rd:
                if hasattr(node, 'rd_map'):
                    rd_label = ete3.faces.TextFace(
                        str(node.rd_map) if hasattr(node, 'rd_map') else (0),
                        fgcolor='Green')
                    if node.rd_map == 100:
                        rd_label.inner_border.type=0
                        rd_label.inner_border.width=1
                    node.add_face(rd_label, column=0)
            if map_iqtree:
                if hasattr(node, 'iq_map'):
                    iq_label = ete3.faces.TextFace(
                        str(node.iq_map) if hasattr(node, 'iq_map') else (0),
                        fgcolor='Red')
                    if node.iq_map == 100:
                        iq_label.inner_border.type=0
                        iq_label.inner_border.width=1
                    node.add_face(iq_label, column=0)
            if node.name:
                node.add_face(ete3.faces.TextFace(node.name), column = 1)

        ts = ete3.TreeStyle()
        cur_col = 0
        if map_rd:
            rd_leg = ete3.faces.TextFace("RootDigger", fgcolor="Green")
            ts.legend.add_face(rd_leg, column=cur_col)
        if map_iqtree:
            iq_leg = ete3.faces.TextFace("IQ-TREE", fgcolor="Red")
            ts.legend.add_face(iq_leg, column=cur_col)
        ts.layout_fn = layout
        ts.show_leaf_name = False
        tree_image_name = "{taxa}_{align}.png".format(
            taxa=result['tree'], align=result['alignment'])
        mapped_tree.render(tree_image_name, tree_style=ts)


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
    parser.add_argument('--atol', type=float, default=1e-4)
    parser.add_argument('--factor', type=float, default=1e7)
    parser.add_argument('--bfgstol', type=float, default=1e-4)
    parser.add_argument('--run-rd', dest='runrd', action='store_true')
    parser.add_argument('--run-iq-tree', dest='runiq', action='store_true')
    parser.add_argument('--no-run-rd', dest='runrd', action='store_false')
    parser.add_argument('--no-run-iq-tree', dest='runiq', action='store_false')
    parser.set_defaults(runrd=True)
    parser.set_defaults(runiq=True)
    args = parser.parse_args()

    RD += " --atol {atol} --factor {factor} --bfgstol {bfgstol}".format(
        atol=args.atol, factor=args.factor, bfgstol=args.bfgstol)

    if args.runiq and not shutil.which("iqtree"):
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
            experiments.append(
                exp('.', i, trees, aligns, args.runrd, args.runiq))

        PROGRESS_BAR.update(PROGRESS_BAR_ITER.value)
        PROGRESS_BAR_ITER.value += 1
        with multiprocessing.Pool(args.procs) as tp:
            finished_exp = tp.map(exp.run_all, experiments)
        experiment_summary = summary(finished_exp)
        experiment_summary.write('test_results')
        produce_mapped_root_images('test_results.json', args.runiq, args.runrd)
