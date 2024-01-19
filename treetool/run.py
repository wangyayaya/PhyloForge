#!/usr/bin/env python
# -*- coding:utf-8 -*-


import argparse
import os.path
import sys
import configparser
import treetool


def get_parser():
    parser = argparse.ArgumentParser(description=f'''
    TreeTool: A tool for multiple phylogenetic analyses
    ***** https://github.com/wangyayaya/TreeTool *****''')

    parser.add_argument("-l", "--lcn", type=str, dest='lcn_tree',
                        help="construct phylogenetic tree based on low-copy nuclear genes ")
    parser.add_argument("-m", "--mul", type=str, dest='mul_tree',
                        help="construct phylogenetic tree based on low-copy nuclear genes (retain multiple copies)")
    parser.add_argument("-s", "--snp", type=str, dest='snp_tree', help="construct phylogenetic tree based on SNP sites")
    parser.add_argument("-S", "--sv", type=str, dest='sv_tree',
                        help="construct phylogenetic tree based on SV (structural variation)")
    parser.add_argument("-o", "--organelle", type=str, dest='organelle_tree',
                        help="construct phylogenetic tree based on organelle-encoded genes.")
    parser.add_argument("-w", "--whole_genome", type=str, dest='whole_genome_tree',
                        help="construct phylogenetic tree based on organelle whole genome sequence.")
    parser.add_argument("-g", "--gene", type=str, dest='gene_tree', help="construct gene tree")
    parser.add_argument("-c", action="store_true",
                        help="get the run configuration file: treetool -c >run.cfg")
    parser.add_argument("-gsc", action="store_true",
                        help="get the software configuration file: treetool -gsc >software.cfg")
    parser.add_argument("-sc", help="perform software configuration: treetool -sc software.cfg")
    parser.add_argument("-v", "--version", action="version", version="treetool 测试版")
    args = parser.parse_args()

    if args.c:
        cfg = os.path.join(treetool.__path__[0], "cfg/run.config")
        with open(cfg) as f:
            for line in f.readlines():
                print(line, end='')
        sys.exit()

    elif args.gsc:
        cfg = os.path.join(treetool.__path__[0], "cfg/software.cfg")
        with open(cfg) as f:
            for line in f.readlines():
                print(line, end='')
        sys.exit()

    elif args.lcn_tree:
        opt_cfg = args.lcn_tree
        tree = 'lcn'

    elif args.mul_tree:
        opt_cfg = args.mul_tree
        tree = 'mul'

    elif args.snp_tree:
        opt_cfg = args.snp_tree
        tree = 'snp'

    elif args.sv_tree:
        opt_cfg = args.sv_tree
        tree = 'sv'

    elif args.organelle_tree:
        opt_cfg = args.organelle_tree
        tree = 'organelle'

    elif args.whole_genome_tree:
        opt_cfg = args.whole_genome_tree
        tree = 'whole_genome'

    elif args.gene_tree:
        opt_cfg = args.gene_tree
        tree = 'gene'

    elif args.sc:
        conf = configparser.ConfigParser()
        conf.read(args.sc)
        cfg = os.path.join(treetool.__path__[0], "cfg/software.ini")
        existing_conf = configparser.ConfigParser()
        existing_conf.read(cfg)
        if conf.has_section('software'):
            for section in conf.sections():
                if not existing_conf.has_section(section):
                    existing_conf.add_section(section)
                for option in conf.options(section):
                    value = conf.get(section, option)
                    existing_conf.set(section, option, value)
            with open(cfg, 'w') as f:
                existing_conf.write(f)
            print('Software config file has been updated')
        else:
            print('Invalid config file format')
        sys.exit()


    else:
        print(parser.format_usage())
        sys.exit()
    # opt作为参数传入其它模块，这样其它模块中就不用引入run？
    return [opt_cfg, tree]


def main():
    tree = get_parser()[1]
    if tree == 'lcn' or tree == 'mul':
        import treetool.lcntree
        treetool.lcntree.LcnTree().lcn_tree()
    elif tree == 'snp':
        import treetool.vcftree
        treetool.vcftree.VcfTree().snp_tree()
    elif tree == 'sv':
        import treetool.svtree
        treetool.svtree.SvTree().sv_tree()
    elif tree == 'organelle':
        import treetool.organelletree
        treetool.organelletree.OrganelleTree().run_organelle_tree()
    elif tree == 'whole_genome':
        import treetool.wholegenometree
        treetool.wholegenometree.WholeGenomeTree().run_whole_genome_tree()
    elif tree == 'gene':
        import treetool.genetree
        treetool.genetree.GeneTree().run_gene_tree()
    else:
        exit()
