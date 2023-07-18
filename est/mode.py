"""
运行步骤
0：从头跑完全流程
1：从头跑到筛选低拷贝OG，初步判断OG数目是否适当
3：只运行筛选OG这一步骤，通过改变cover及copy_number参数改变OG数量，直到满足自己要求
4：OG数量恰当后，即跑完步骤3后，开始跑完后续全部步骤
"""
import sys

'''from est import script, hmm2OG'''
import script, hmm2OG
cmd = script.RunCmd()


class MODE:
    def __init__(self):
        self.seq = 0
        self.coa_con = 0
        self.mode = 0
        opt = hmm2OG.HMM_OG().get_config()
        for k, v in opt.items():
            setattr(self, str(k), v)

    def run_mode0(self):
        """cds to species tree"""
        statuss= cmd.run_hmmscan_pl()
        # statuss = [0]
        if 1 in statuss:
            print("The hmmsearch program failed to run")
            sys.exit(1)
        else:
            hmm2OG.HMM_OG().run_hmm2OG()
            cmd.run_genetree_mul()
            cmd.run_astral()
            if self.coa_con == 1:
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
        """select OG"""
        hmm2OG.HMM_OG().run_hmm2OG()

    def run_mode3(self):
        """OG to tree"""
        cmd.run_genetree_mul()
        cmd.run_astral()
        if int(self.coa_con) == 1:
            cmd.run_contree()

    def run_mode4(self):
        """重新选择做树软件"""
        pass

    def run(self):
        try:
            n = int(self.mode)
            if n == 0:
                run = self.run_mode0()
            elif n ==1:
                run = self.run_mode1()
            elif n ==2:
                run = self.run_mode2()
            elif n ==3:
                run = self.run_mode3()
            elif n == 4:
                run =self.run_mode4()
            else:
                print("please check")

        except ValueError:
            print("please check type of mode (int)")
            sys.exit()

        return run
