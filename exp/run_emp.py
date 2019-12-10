#!/usr/bin/env python3

import multiprocessing.pool
import subprocess
import yaml
import os
import datetime

RD = os.path.abspath(
    "../bin/rd"
) + " --msa {msa} --tree {tree} --exhaustive --early-stop --treefile {treefile}"


def produce_rd_command(msa, tree, treefile):
    return RD.format(msa=msa, tree=tree, treefile=treefile).split(' ')


def make_done_file(msa):
    return os.path.join(os.path.dirname(msa),
                        "." + os.path.basename(msa) + ".done")


def clean_done(msa):
    os.remove(make_done_file(msa))

def set_done(msa):
    with open(make_done_file(msa), 'w') as df:
        df.write(datetime.datetime.now().isoformat())


def check_done(msa):
    return os.path.exists(make_done_file(msa))


def run_rd(msa, tree, treefile, log):
    with open(log, 'a') as logfile:
        subprocess.run(produce_rd_command(msa, tree, treefile), stdout=logfile)
    set_done(msa)


def make_file_paths(directory_root, dataset):
    msa = os.path.join(directory_root, ds['msa'])
    tree = os.path.join(directory_root, ds['tree'])
    results = os.path.join(directory_root, ds['results'])
    log = os.path.join(directory_root, ds['log'])
    return (msa, tree, results, log)


def verify_dataset(directory_root, datasets):
    for ds in datasets:
        msa, tree, results, log = make_file_paths(directory_root, ds)
        if not os.path.exists(msa):
            raise Exception("msa file: " + msa + " does not exist")
        if not os.path.exists(tree):
            raise Exception("tree file: " + tree + " does not exist")


with open('./datasets.yaml') as infile:
    datasets = yaml.load(infile, Loader=yaml.FullLoader)

jobs = []

for k, d in datasets.items():
    directory_root = d['directory']
    for ds in d['datasets']:
        verify_dataset(directory_root, ds)
        jobs.append(make_file_paths(directory_root, ds))

tp = multiprocessing.pool.ThreadPool()
tp.starmap(run_rd, jobs)
