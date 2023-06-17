import argparse
import os
import sys


from script import script,hmm2OG,mode
cmd = script.RunCmd()


if __name__ == '__main__':
    #cmd.mkdir()
    #script.format_and_trans_mul()
    #cmd.run_hmmscan_pl()
    #hmm2OG.HMM_OG().run_hmm2OG()
    #cmd.run_tree_mul()
    #cmd.get_super_gene()
    mode.MODE().mode0()
