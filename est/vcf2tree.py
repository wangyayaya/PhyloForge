"""snp vcf to tree"""
import os
import subprocess
import sys
from configparser import ConfigParser


class VCF2TREE():
    def __init__(self):
        self.treesoft = 'treebest'

    def run_command(self, command):
        try:
            output = subprocess.check_output(command, shell=True, stderr=subprocess.DEVNULL)
            return 0, output.decode()
        except subprocess.CalledProcessError as e:
            print(e)
            return e.returncode

    def convert_vcf_to_seq(self, vcf_file, out, f):
        # Open the input vcf file.
        with open(vcf_file) as vcf_file:
            lines = vcf_file.readlines()

        # Extract sequence data from vcf file.
        ids = []
        seqs = []
        ids_found = False
        for line in lines:
            if line.startswith("#CHROM"):
                ids = line.strip().split()[9:]
                seqs = [""] * len(ids)
                ids_found = True
            elif ids_found and line.strip():
                fields = line.strip().split()
                ref = fields[3]
                alleles = [ref] + fields[4].split(",")
                for i, gt in enumerate(fields[9:], start=0):
                    if gt.split(":")[0] == "./." or gt.split(":")[0] == ".|.":
                        gt_seq = "-"
                    else:
                        if '/' in gt:
                            alleles_indices = [int(idx) for idx in gt.split("/") if idx.isdigit()]
                        elif '|' in gt:
                            alleles_indices = [int(idx) for idx in gt.split("|") if idx.isdigit()]
                        gt_seq = "".join(alleles[idx] for idx in alleles_indices if idx < len(alleles))
                        for nu in gt_seq:
                            if nu not in ['A', 'T', 'G', 'C', 'a', 't', 'g', 'c', '-']:
                                gt_seq = gt_seq.replace(nu, '-')
                    # VCF4.2格式中，同一位点突变可能会涉及多个碱基，以最长的那一个样本的碱基长度作为所有样本长度，不足的用'-'补齐
                    len_gt = 1
                    for allele in alleles:
                        if len(allele) > len_gt:
                            len_gt = len(allele)
                    sup_lenth = len_gt - len(gt_seq)
                    gt_seq = gt_seq + '-' * sup_lenth

                    seqs[i] += gt_seq

        if f == 'phy':
            # Prepare the phylip string.
            phylip_str = f"{len(ids)} {len(seqs[0])}\n"
            for i in range(len(ids)):
                phylip_str += f"{ids[i].ljust(10)} {seqs[i]}\n"
            out_str = phylip_str
        # 60个字符换一行
        else:
            fasta_str = ""
            for i in range(len(ids)):
                fasta_str += f">{ids[i]}\n"
                seq = seqs[i]
                while len(seq) > 0:
                    fasta_str += seq[:60] + "\n"
                    seq = seq[60:]  # 每一次删掉序列前60个碱基
            out_str = fasta_str
        # fa格式不换行
        """else:
            fasta_str = ""
            for i in range(len(ids)):
                fasta_str += f">{ids[i]}\n{seqs[i]}\n"
            out_str = fasta_str"""

        # Write the phylip string to file.
        with open(out, "w") as out:
            out.write(out_str)

    def PHYLIP_tree(self, infile, out_path):
        # 这个软件在运行过程中会一直在屏幕打印信息，用什么系统运行命令都会有问题，暂时不考虑用这个软件
        seqboot = f'{infile}\nr\n100\ny\n9'
        dnadist = 'seqboot.out\nt\n2.3628\nm\nd\n100\n2\ny'
        neighbor = 'dnadist.out\nm\n100\n9\ny'
        consense = 'nei.tree\ny'

        seqboot_path = out_path + '/seqboot.par'
        seqboot_par = open(seqboot_path, 'w')

        dnadist_path = out_path + '/dnadist.par'
        dnadist_par = open(dnadist_path, 'w')

        neighbor_path = out_path + '/neighbor.par'
        neighbor_par = open(neighbor_path, 'w')

        consense_path = out_path + '/consense.par'
        consense_par = open(consense_path, 'w')

        print(seqboot, file=seqboot_par)
        print(dnadist, file=dnadist_par)
        print(neighbor, file=neighbor_par)
        print(consense, file=consense_par)

        run = f'cd {out_path} && seqboot< seqboot.par && mv outfile seqboot.out &&dnadist<dnadist.par && ' \
              'mv outfile dnadist.out && neighbor<neighbor.par && mv outfile nei.out && mv outtree nei.tree && ' \
              'consense<consense.par && mv outfile cons.out && mv outtree constree'
        return self.run_command(run)

    def PhyML_tree(self, in_phy):
        run = f'{self.treesoft} -i {in_phy} -b 100 -m HKY85 -f m -v e -a e -o tlr'
        return self.run_command(run)

    def tree_best_tree(self, in_fa):
        run = f'{self.treesoft} nj -b 1000  in_fa >{in_fa}_treebest.NHX'
        return self.run_command(run)


def main():


if __name__ == '__main__':
    VCF2TREE().tree_best_tree('a')
