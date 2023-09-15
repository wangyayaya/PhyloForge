#!/usr/bin/env python
# -*- coding:utf-8 -*-


import argparse
import os.path
import sys
import configparser

import est


def get_parser():
    parser = argparse.ArgumentParser(description=f'''
    EST (Easy species tree): a user-friendly pipeline used to construct phylogenetic tree based on low-copy nucleic acid (LCN) or SNP sites.
    ***** https://github.com/wangyayaya/EasySpeciesTree *****''')

    parser.add_argument("-l", "--lcn", type=str, dest='lcn_tree', help="construct phylogenetic tree based on low copy "
                                                                       "genes, you can obtain configfile by running: "
                                                                       "est -grc >run.cfg")
    parser.add_argument("-s", "--snp", type=str, dest='snp_tree', help="construct phylogenetic tree based on SNP sites")
    parser.add_argument("-o", "--organelle", type=str, dest='organelle_tree',
                        help="construct phylogenetic tree based on organelle-encoded genes.")
    parser.add_argument("-grc", action="store_true", help="get the run configuration file")
    parser.add_argument("-gsc", action="store_true", help="get the software configuration file")
    parser.add_argument("-sc", help="perform software configuration, you can obtain configfile by running: est -gsc "
                                    ">software.cfg")
    parser.add_argument("-v", "--version", action="version", version="est 测试版")
    args = parser.parse_args()

    if args.grc:
        cfg = os.path.join(est.__path__[0], "cfg/run.config")
        with open(cfg) as f:
            for line in f.readlines():
                print(line, end='')
        sys.exit()

    elif args.gsc:
        cfg = os.path.join(est.__path__[0], "cfg/software.cfg")
        with open(cfg) as f:
            for line in f.readlines():
                print(line, end='')
        sys.exit()

    elif args.lcn_tree:
        opt_cfg = args.lcn_tree
        tree = 'lcn'

    elif args.snp_tree:
        opt_cfg = args.snp_tree
        tree = 'snp'

    elif args.organelle_tree:
        opt_cfg = args.organelle_tree
        tree = 'organelle'

    elif args.sc:
        conf = configparser.ConfigParser()
        conf.read(args.sc)
        cfg = os.path.join(est.__path__[0], "cfg/software.ini")
        if conf.sections()[0] == 'software':
            conf.write(open(cfg, 'w'))
            print('Software config file has been modified')
        else:
            print('Software config file no change')
        sys.exit()

    else:
        print(parser.format_usage())
        sys.exit()
    # opt作为参数传入其它模块，这样其它模块中就不用引入run这个包？
    return [opt_cfg, tree]


def main():
    import est.lcntree
    import est.vcftree
    import est.organelletree
    tree = get_parser()[1]
    if tree == 'lcn':
        est.lcntree.LcnTree().lcn_tree()
    elif tree == 'snp':
        est.vcftree.VcfTree().snp_tree()
    elif tree == 'organelle':
        est.organelletree.OrganelleTree().run_organelle_tree()
    else:
        sys.exit()
