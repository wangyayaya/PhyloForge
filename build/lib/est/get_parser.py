import argparse
import os
import sys


def get_parser():
    parser = argparse.ArgumentParser(description=f"To get started: {__file__} -c run.config")
    parser.add_argument("-c", "--config", type=str, help=f"Run main program, You can obtain configfile by running: "
                                                         f"{__file__} -g >run.cfg, or refer to this file: "
                                                         f"{os.path.splitext(os.path.split(__file__)[0])[0]}/run.config")
    parser.add_argument("-g", "--getcfg", action="store_true", help="Get the config file")
    args = parser.parse_args()
    if args.getcfg:
        print("""
[opt]
in_path = /path/to/cds_dir/
out_path = /path/to/output/
thread = 1
orthodb = /path/to/buscoDB
aln_software = muscle
tree_software = iqtree
cover = 7
copy_number = 1
step = 0
seq = cds
coa_con = 0

[software]
hmmscan=/path/to/hmmscan
mafft=/path/to/mafft
muscle=/path/to/muscle
trimal=/path/to/trimal
raxmal=/path/to/raxmal
iqtree=/path/to/iqtree
astral=/path/to/astral.jar     
        """)
        sys.exit()

    elif args.config:
        opt_cfg = args.config

    else:
        print(parser.format_usage())
        sys.exit()

    return opt_cfg

