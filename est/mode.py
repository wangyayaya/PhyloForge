"""
运行步骤
0：从头跑完全流程
1：从头跑到筛选低拷贝OG，初步判断OG数目是否适当
3：只运行筛选OG这一步骤，通过改变cover及copy_number参数改变OG数量，直到满足自己要求
4：OG数量恰当后，即跑完步骤3后，开始跑完后续全部步骤
5：重新构树，即重新选择序列比对即构树软件
"""
import sys

from est import script, hmm2OG

cmd = script.RunCmd()


class MODE:
    def __int__(self):
        self.seq = 0
        self.coa_con = 0
        opt = hmm2OG.HMM_OG().get_config()
        for k, v in opt.items():
            setattr(self, str(k), v)

    def mode0(self):
        """从头开始运行全流程cds to species tree"""
        # statuss = cmd.run_hmmscan_pl()
        statuss = [0]
        if 1 in statuss:
            print("The hmmsearch program failed to run")
            sys.exit(1)
        else:
            hmm2OG.HMM_OG().run_hmm2OG()
            statuss_aln, statuss_trim, statuss_built_tree = cmd.run_genetree_mul()
            if any([status == 1 for status in statuss_aln + statuss_trim + statuss_built_tree]):
                print("hhh")
            else:
                cmd.run_astral()

    def mode1(self):
        """从头运行至筛选OG,cds to SOG/LOG"""
        cmd.mkdir()
        cmd.run_format_and_trans()
        statuss = cmd.run_hmmscan_pl()
        if 1 in statuss:
            print("The hmmsearch program failed to run")
            sys.exit(1)
        else:
            hmm2OG.HMM_OG().run_hmm2OG()

    def mode2(self):
        """只运行筛选OG这一步骤"""
        hmm2OG.HMM_OG().run_hmm2OG()

    def mode3(self):
        """筛选出适当的OG数目后，跑完后续全部流程"""
        cmd.run_genetree_mul()
        cmd.run_astral()

    def mode4(self):
        """重新选择做树软件"""
        pass


def run():
    pass
