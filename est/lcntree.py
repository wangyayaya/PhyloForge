"""
mode to run
0：from cds to tree
1：from cds to OGs
2：select OGs
3：from OGs to tree
"""
import sys

from configparser import ConfigParser

from est import script, hmm2OG, run

# import script, hmm2OG
cmd = script.RunCmd()


class LcnTree:
    def __init__(self):
        self.coa_con = 0
        self.mode = 0
        opt_cfg = run.get_parser()[0]
        opt = cmd.get_config(opt_cfg, 'lcn_opt')
        #opt = self.get_config('lcn_opt')
        for k, v in opt.items():
            setattr(self, str(k), v)

    def check(self):
        global n
        try:
            int(self.cover)
        except AttributeError:
            pass
        except ValueError:
            print('Please check the cover, which must be a int')
            sys.exit()

        try:
            int(self.copy_number)
        except AttributeError:
            pass
        except ValueError:
            print('Please check the copy_number, which must be a int')
            sys.exit()

        try:
            if int(self.coa_con) not in [0, 1]:
                print('coa_con should be 0 or 1, please check it')
                sys.exit()
        except AttributeError:
            pass
        except ValueError:
            print("coa_con should be 0 or 1, please check it")
            sys.exit()

        try:
            if str(self.seq).upper() not in ["CDS", "PEP", "CODON"]:
                print(f'The seq: {self.seq} is not recognized, and CDS will be used')
        except AttributeError:
            print('The seq is not specified, and CDS will be used')

        try:
            int(self.thread)
        except AttributeError:
            pass
        except ValueError:
            print("thread should be int (>=2), please check it")
            sys.exit()

        try:
            if str(self.aln_software).upper() not in ['MAFFT', 'MUSCLE']:
                print(f'The aln_software: {self.aln_software} is not recognized, and MAFFT will be used or sequence alignment')
        except AttributeError:
            print('Sequence alignment software is not specified, muscle software will be used for sequence alignment')

        try:
            if str(self.tree_software).upper() not in ['RAXML', 'IQTREE', 'FASTTREE']:
                print(f'The tree_software: {self.tree_software} is not recognized, and RAxML will be used')
        except AttributeError:
            print('The tree_software is not specified, raxml software will be used')

        try:
            n = int(self.mode)
            if  n not in [0, 1, 2, 3]:
                print('mode should be 0, 1, 2 or 3, please check it')
                sys.exit()
        except ValueError:
            print("please check type of mode (int)")
            sys.exit()
        except AttributeError:
            pass

        return n


    def run_mode0(self):
        """cds to species tree"""
        statuss = cmd.run_hmmscan_pl()
        # statuss = [0]
        if 1 in statuss:
            print("The hmmsearch program failed to run")
            sys.exit(1)
        else:
            hmm2OG.HMM_OG().run_hmm2OG()
            if int(self.coa_con) == 0:
                cmd.run_genetree_mul()
                cmd.run_astral()
            elif int(self.coa_con) == 1:
                cmd.run_genetree_mul()
                cmd.run_astral()
                cmd.run_contree()
            else:
                print("The parameter coa_con must be selected between 1 and 2!")
                sys.exit()

    def run_mode1(self):
        """cds to SOG/LOG"""
        statuss = cmd.run_hmmscan_pl()
        if 1 in statuss:
            print("The hmmsearch program failed to run")
            sys.exit(1)
        else:
            hmm2OG.HMM_OG().run_hmm2OG()

    def run_mode2(self):
        hmm2OG.HMM_OG().run_hmm2OG()

    def run_mode3(self):
        """OG to tree"""
        if int(self.coa_con) == 0:
            cmd.run_genetree_mul()
            cmd.run_astral()
        # elif int(self.coa_con) == 1:
        #    cmd.run_contree()
        elif int(self.coa_con) == 1:
            cmd.run_genetree_mul()
            cmd.run_astral()
            cmd.run_contree()

    def lcn_tree(self):
        n = self.check()
        if n == 0:
            self.run_mode0()
        elif n == 1:
            self.run_mode1()
        elif n == 2:
            self.run_mode2()
        elif n == 3:
            self.run_mode3()
        else:
            self.run_mode0()
