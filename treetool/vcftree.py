"""snp vcf to tree"""
import os
from treetool import script
from treetool.script import RunCmd
cmd = script.RunCmd()


class VcfTree(RunCmd):
    def __init__(self):
        super().__init__()

    def generate_ambiguous_code(self, base1, base2):
        # 非纯和突变碱基的歧义码IUPAC Ambiguity Codes，类型没有vcf2phylip的多
        # 但运行结果检查了一下和那个脚本是一样的，因为只考虑单碱基突变的话就只有下面这几种突变
        # 一个自己的想法，如果一个位点参考是A，突变是T，那么存在杂合突变时是否可以将0/0转换为AA，1/1转换为TT，而0/1转换为AT，感觉可行
        base1 = base1.upper()
        base2 = base2.upper()
        # 这里是获得的基因型，ATGC*.,"."是缺失的情况
        # 不考虑一个位点涉及多个碱基的情况
        if base1 == base2:
            if {base1} <= {"A", "T", "G", "C"}:
                return base1
            elif base1 == '.':
                return "-"
            elif base1 == '*':
                return 'N'
        elif base1 not in ["A", "T", "G", "C"] and base2 not in ["A", "T", "G", "C"]:
            return "N"
        elif {base1, base2} <= {"A", "C"}:
            return "M"
        elif {base1, base2} <= {"A", "G"}:
            return "R"
        elif {base1, base2} <= {"A", "T"}:
            return "W"
        elif {base1, base2} <= {"C", "G"}:
            return "S"
        elif {base1, base2} <= {"C", "T"}:
            return "Y"
        elif {base1, base2} <= {"G", "T"}:
            return "K"
        elif {base1, base2} <= {"A", "*"}:
            return "a"
        elif {base1, base2} <= {"T", "*"}:
            return "t"
        elif {base1, base2} <= {"G", "*"}:
            return "g"
        elif {base1, base2} <= {"C", "*"}:
            return "c"
        else:
            return "N"

    def convert_vcf_to_seq(self, vcf_file, out_file, format, LF):
        if LF.upper().startswith("F"):
            LF = False
        else:
            LF = True
        # Open the input vcf file.
        with open(vcf_file) as vcf_file:
            lines = vcf_file.readlines()
        for line in lines:
            if line.startswith("#CHROM"):
                sample_name = line.strip().split()[9:]
                # print(sample_name)
                seqs = [""] * len(sample_name)
                # print(seqs)
            elif line.startswith("#"):
                pass
            elif line.strip():
                fields = line.strip().split()
                ref = fields[3]
                alleles = [ref] + fields[4].split(",")
                # print(alleles)
                # 这里只考虑每个位点只涉及单个碱基突变的情况，多个碱基突变的话就很难考虑那个歧义码。之前写的脚本是以最长那一个突变为准
                # 不考虑歧义码，其它缺失的用‘-’来补齐，如果同时考虑歧义码及多个碱基，没办法转换
                base_len_le_1 = True
                for base in alleles:
                    base_len = len(base)
                    if base_len > 1:
                        base_len_le_1 = False
                        break
                if base_len_le_1:
                    for i, alleles_indices in enumerate(fields[9:], start=0):
                        # i是转化后sample_name列表中物种名的索引
                        if alleles_indices.split(":")[0]:
                            alleles_indices = alleles_indices.split(":")[0]
                            try:
                                if '/' in alleles_indices:
                                    allele1 = alleles_indices.split('/')[0]
                                    allele2 = alleles_indices.split('/')[1]
                                elif '|' in alleles_indices:
                                    allele1 = alleles_indices.split('|')[0]
                                    allele2 = alleles_indices.split('|')[1]
                                base1 = alleles[int(allele1)]
                                base2 = alleles[int(allele2)]
                            except ValueError:
                                base1 = base2 = '.'
                            # print(base1, base2)
                            base = self.generate_ambiguous_code(base1, base2)
                            # print(base)
                            seqs[i] += base
        # for i in range(len(sample_name)):
        #   print(f">{sample_name[i]}\n{seqs[i]}")

        if format == 'phylip':
            # phy序列不换行
            phylip_str = f"{len(sample_name)} {len(seqs[0])}\n"
            for i in range(len(sample_name)):
                phylip_str += f"{sample_name[i].ljust(10)} {seqs[i]}\n"
            out_str = phylip_str

        elif format == 'fasta' and LF:
            # fa格式60个字符换一行
            fasta_str_LF = ""
            for i in range(len(sample_name)):
                fasta_str_LF += f">{sample_name[i]}\n"
                seq = seqs[i]
                while len(seq) > 0:
                    fasta_str_LF += seq[:60] + "\n"
                    seq = seq[60:]  # 每一次删掉序列前60个碱基
            out_str = fasta_str_LF

        elif format == 'fasta' and not LF:
            # fa格式不换行
            fasta_str = ""
            for i in range(len(sample_name)):
                fasta_str += f">{sample_name[i]}\n{seqs[i]}\n"
            out_str = fasta_str

        with open(out_file, 'w') as f:
            f.writelines(out_str)

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
        return cmd.run_command(run)

    def snp_tree(self):
        basename = 'SNP'
        if self.tree_software.upper() == 'TREEBEST':
            self.convert_vcf_to_seq(self.vcf_file, f'{self.out_path}/SNP.fa', 'fasta', 'False')
            infile = f'{self.out_path}/SNP.fa'
        else:
            self.convert_vcf_to_seq(self.vcf_file, f'{self.out_path}/SNP.phy', 'phylip', 'False')
            infile = f'{self.out_path}/SNP.phy'
        cmd.built_tree(infile, f'{self.out_path}', basename, 1)


