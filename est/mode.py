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
opt_cfg = run.get_parser()[0]


class MODE:
    def __init__(self):
        self.seq = 0
        self.coa_con = 0
        self.mode = 0
        opt = self.get_config()
        for k, v in opt.items():
            setattr(self, str(k), v)

    def get_config(self):
        cfg_parser = ConfigParser()
        # read options
        cfg_parser.read(opt_cfg)
        opt = dict(cfg_parser.items('lcn_opt'))
        return opt

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
                cmd.run_contree()
            else:
                cmd.run_genetree_mul()
                cmd.run_astral()
                cmd.run_contree()

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
        elif int(self.coa_con) == 1:
            cmd.run_contree()
        elif int(self.coa_con) == 2:
            cmd.run_genetree_mul()
            cmd.run_astral()
            cmd.run_contree()

    def run(self):
        try:
            n = int(self.mode)
        except ValueError:
            print("please check type of mode (int)")
            sys.exit()

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


