#!/usr/bin/env python3

import ete3
import math
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--tree', type=str)
parser.add_argument('--image', type=str)
parser.add_argument('--draw-names', action='store_true')
parser.add_argument('--root', type=str)

color1 = "#1f77b4"
color2 = "#ff7f0e"
branch_color = "#7f7f7f"
root_color = '#d62728'

def make_render_node(draw_names):
    def render_node(node):
        cur_col = 0
        ns = ete3.NodeStyle()
        ns['fgcolor'] = branch_color
        if hasattr(node, 'name') and draw_names:
            label = ete3.faces.TextFace(node.name, fgcolor='black')
            node.add_face(label, column=cur_col)
            cur_col += 1
        if hasattr(node, 'LWR'):
            lwr = float(node.LWR)
            label = ete3.faces.TextFace("{:.4}".format(lwr), fgcolor=color1)
            if not abs(lwr) < 1e-7:
                node.add_face(label, column=cur_col)
        if hasattr(node, 'alpha'):
            alpha = float(node.alpha)
            label = ete3.faces.TextFace("{:.3}".format(alpha), fgcolor=color2)
            if not abs(lwr) < 1e-7:
                node.add_face(label, column=cur_col)
        if hasattr(node, 'root') and node.root:
            ns['fgcolor'] = root_color
            ns['hz_line_color'] = root_color
            ns['hz_line_width'] = 4
        if hasattr(node, 'foo'):
            label = ete3.faces.TextFace(node.foo, fgcolor='Green')
            node.add_face(label, column=cur_col)
        if hasattr(node, 'fizz'):
            label = ete3.faces.TextFace(node.fizz, fgcolor='Red')
            node.add_face(label, column=cur_col)
        node.set_style(ns)
    return render_node


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

def mark_root(tree, rooted_tree):
    extract_node_with_clade_pair(tree, get_root_clades(rooted_tree)).root = True

def get_leaf_set(tree):
    ls = []
    for n in tree.traverse():
        if n.is_leaf():
            ls.append(str(n))
    return sorted(ls)

def draw_lh(root, image, tree, draw_names):
    ts = ete3.TreeStyle()

    if not os.path.isfile(tree):
        raise Exception("File not found:" + tree)

    with open(tree) as infile:
        t = ete3.Tree(infile.read())

    if root:
        with open(root) as root_tree_file:
            root_tree = ete3.Tree(root_tree_file.read())
            if get_leaf_set(t) != get_leaf_set(root_tree):
                raise Exception("leaves don't match")
            mark_root(t, root_tree)


    ts.layout_fn = make_render_node(draw_names)
    leg1 = ete3.faces.TextFace("LWR", fgcolor=color1)
    leg2 = ete3.faces.TextFace("Alpha", fgcolor=color2)
    leg3 = ete3.faces.TextFace("Root Branch", fgcolor=root_color)
    ts.legend.add_face(leg1, column=0)
    ts.legend.add_face(leg2, column=0)
    if root:
        ts.legend.add_face(leg3, column=0)
    ts.show_leaf_name = False

    t.render(image, tree_style=ts)

if __name__ == "__main__":
    args = parser.parse_args()
    draw_names = args.draw_names
    draw_lh(args.root, args.image, args.tree)
