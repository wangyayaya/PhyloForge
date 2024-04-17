"""
mode to run
0：from cds to tree
1：from cds to OGs
2：select OGs
3：from OGs to tree
"""
import sys
import os
import datetime
from treetool import script, hmm2OG
from treetool.script import RunCmd

cmd = script.RunCmd()


class LcnTree(RunCmd):
    def __init__(self):
        self.mode = 0
        super().__init__()
        if self.retain_multi_copy == 'True':
            self.retain_multi_copy = True
            self.coa_con = 0
            try:
                os.path.abspath(os.path.dirname(self.astral_pro))
            except AttributeError:
                current_time = datetime.datetime.now()
                print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] If [retain_multi_copy = Ture], astral-pro "
                      f"software is not configured, please refer to the README file for configuration.")
                sys.exit(1)
        else:
            self.retain_multi_copy = False

    def HMM_OG(self):
        if self.retain_multi_copy:
            hmm2OG.HMM_OG().run_hmm2OG_multi_copy()
        else:
            hmm2OG.HMM_OG().run_hmm2OG()

    def ortho_OG(self):
        if not self.retain_multi_copy:
            hmm2OG.HMM_OG().ortho2OG(False)
        else:
            hmm2OG.HMM_OG.ortho2OG(True)

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
            n = int(self.mode)
            if n not in [0, 1, 2, 3]:
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
        cmd.run_hmmscan_pl()
        self.HMM_OG()
        if int(self.coa_con) == 0:
            cmd.run_genetree_mul()
            cmd.run_astral()
        elif int(self.coa_con) == 1:
            cmd.run_genetree_mul()
            cmd.run_astral()
            if not self.retain_multi_copy:
                cmd.run_contree()
        else:
            print("The parameter coa_con must be selected between 0 and 1!")
            sys.exit()

    def run_mode1(self):
        """cds to SOG/LOG"""
        cmd.run_hmmscan_pl()
        self.HMM_OG()

    def run_mode2(self):
        self.HMM_OG()

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
            sys.exit()
