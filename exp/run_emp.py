#!/usr/bin/env python3

import multiprocessing.pool
import subprocess
import yaml
import os
import datetime
import importlib
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--prefix', type=str, required=True)
args = parser.parse_args()

drawlh = importlib.import_module('drawlh')

RD = os.path.abspath(
    "../bin/rd"
) + " --msa {msa} --tree {tree} --exhaustive --treefile {treefile} {options}"

def produce_rd_command(msa, tree, treefile, early_stop):
    options = ''
    if early_stop:
        options += '--early-stop'
    return RD.format(msa=msa, tree=tree, treefile=treefile, options=options).split(' ')


def make_done_file(log):
    return os.path.join(os.path.dirname(log),
                        "." + os.path.basename(log) + ".done")


def clean_done(log):
    os.remove(make_done_file(log))

def set_done(log):
    with open(make_done_file(log), 'w') as df:
        df.write(datetime.datetime.now().isoformat())


def check_done(log):
    return os.path.exists(make_done_file(log))


def run_rd(msa, tree, treefile, image, log, options):
    if check_done(log):
        return
    print("rd: ", log)
    with open(log, 'a') as logfile:
        subprocess.run(produce_rd_command(msa, tree, treefile,
            options['early-stop']), stdout=logfile)
    set_done(log)


def make_file_paths(directory_root, prefix, dataset):
    msa = os.path.join(directory_root, ds['msa'])
    tree = os.path.join(directory_root, ds['tree'])
    results = os.path.join(prefix, directory_root, ds['results'])
    image = os.path.join(prefix, directory_root, ds['image'])
    log = os.path.join(prefix, directory_root, ds['log'])
    results_path = os.path.join(prefix, directory_root)
    if not os.path.exists(results_path):
        os.mkdir(results_path)
    return (msa, tree, results, image, log, ds['options'])


def verify_dataset(directory_root,prefix, datasets):
    for ds in datasets:
        msa, tree, results, image, log, early_stop = make_file_paths(directory_root, prefix, ds)
        if not os.path.exists(msa):
            raise Exception("msa file: " + msa + " does not exist")
        if not os.path.exists(tree):
            raise Exception("tree file: " + tree + " does not exist")


with open('./datasets.yaml') as infile:
    datasets = yaml.load(infile, Loader=yaml.FullLoader)

jobs = []
if not os.path.exists(args.prefix):
    os.mkdir(args.prefix)
else:
    if not os.path.isdir(args.prefix):
        raise Exception("Prefix is already a file")


for k, d in datasets.items():
    #directory_root = os.path.join(args.prefix,d['directory'])
    directory_root = d['directory']
    for ds in d['datasets']:
        verify_dataset(directory_root, args.prefix, ds)
        jobs.append(make_file_paths(directory_root,args.prefix, ds))

tp = multiprocessing.pool.ThreadPool()
tp.starmap(run_rd, jobs)

for k, d in datasets.items():
    directory_root = d['directory']
    for ds in d['datasets']:
        try:
            _, tree, results, image, _, _ = make_file_paths(directory_root,
                    args.prefix, ds)
            drawlh.draw_lh(tree, image, results, False)
        except Exception as e:
            print(directory_root,":", ds['msa'], ":", e)
