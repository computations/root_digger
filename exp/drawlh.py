#!/usr/bin/env python3

import ete3
import math
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--tree', type=str)
parser.add_argument('--image', type=str)
args = parser.parse_args()

if not os.path.isfile(args.tree):
    print("file not found")
    os.exit(1)

with open(args.tree) as infile:
    t = ete3.Tree(infile.read())


def render_node(node):
    cur_col = 0
    if hasattr(node, 'name'):
        label = ete3.faces.TextFace(node.name, fgcolor='black')
        node.add_face(label, column=cur_col)
        cur_col += 1
    if hasattr(node, 'LWR'):
        lwr = float(node.LWR)
        label = ete3.faces.TextFace("{:.4}".format(lwr), fgcolor='Green')
        if abs(lwr - 1.0) < 1e-12:
            label.inner_border.type = 0
            label.inner_border.width = 1
        node.add_face(label, column=cur_col)
    if hasattr(node, 'alpha'):
        alpha = float(node.alpha)
        label = ete3.faces.TextFace("{:.3}".format(alpha), fgcolor='Red')
        node.add_face(label, column=cur_col)
    if hasattr(node, 'foo'):
        label = ete3.faces.TextFace(node.foo, fgcolor='Green')
        node.add_face(label, column=cur_col)
    if hasattr(node, 'fizz'):
        label = ete3.faces.TextFace(node.fizz, fgcolor='Red')
        node.add_face(label, column=cur_col)


ts = ete3.TreeStyle()
ts.layout_fn = render_node
leg1 = ete3.faces.TextFace("Likelihood Weight Ratio", fgcolor="Green")
leg2 = ete3.faces.TextFace("Alpha", fgcolor="Red")
ts.legend.add_face(leg1, column=0)
ts.legend.add_face(leg2, column=0)
ts.show_leaf_name = False

t.render(args.image, tree_style=ts)
