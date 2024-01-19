"""snp vcf to tree"""

import sys

from treetool import script
from treetool.script import RunCmd
cmd = script.RunCmd()


class SvTree(RunCmd):
    def __init__(self):
        super().__init__()

    def convert_vcf_to_seq(self, vcf_file, out_file):
        """只对插入和缺失进行转换，插入时：
        如果存在杂合突变：
        ref=0, alt=1, [0,1]
        0/0: 00
        0/1: 01
        1/1: 11
        如果没有杂合突变：
        0/0: 0
        1/1: 1
        缺失时：
        ref=1, alt=0, [1,0]
        如果存在杂合突变：
        0/0: 11
        0/1: 10
        1/1: 00
        如果没有杂合突变：
        0/0: 1
        1/1: 0
        """

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
                # sv_type = fields[4]  符号表示法（Symbolic notation）中第四列即为SV类型
                sv_info = fields[7]
                pairs = sv_info.split(';')
                sv_type = None
                for pair in pairs:
                    key, value = pair.split('=')
                    if key == 'SVTYPE':
                        sv_type = value
                        break
                if sv_type is None or sv_type not in ['DEL', 'INS']:
                    continue
                if sv_type == 'DEL':
                    alleles = [1, 0]
                elif sv_type == 'INS':
                    alleles = [0, 1]
                # print(alleles)
                # 这里好恶心
                # 虽然很恶心，但是改了又报错或者转换的结果不对
                # 恶心就恶心吧，不改了
                Hom = True
                for i, alleles_indices in enumerate(fields[9:], start=0):
                    if alleles_indices.split(":")[0]:
                        alleles_indices = alleles_indices.split(":")[0]
                        # print(alleles_indices)

                        if '/' in alleles_indices:
                            allele1 = alleles_indices.split('/')[0]
                            allele2 = alleles_indices.split('/')[1]
                        elif '|' in alleles_indices:
                            allele1 = alleles_indices.split('|')[0]
                            allele2 = alleles_indices.split('|')[1]

                        if allele1 != allele2:
                            # 只要有一个样本在该位点不是纯和突变，该位点即不是纯和
                            Hom = False
                            break

                if Hom:
                    for i, alleles_indices in enumerate(fields[9:], start=0):
                        if alleles_indices.split(":")[0]:
                            alleles_indices = alleles_indices.split(":")[0]
                            # print(alleles_indices)
                            try:
                                if '/' in alleles_indices:
                                    allele1 = alleles_indices.split('/')[0]
                                    allele2 = alleles_indices.split('/')[1]
                                elif '|' in alleles_indices:
                                    allele1 = alleles_indices.split('|')[0]
                                    allele2 = alleles_indices.split('|')[1]
                                base = alleles[int(allele1)]
                            except ValueError:
                                base = '.'
                        seqs[i] += str(base)
                else:
                    for i, alleles_indices in enumerate(fields[9:], start=0):
                        if alleles_indices.split(":")[0]:
                            alleles_indices = alleles_indices.split(":")[0]
                            # print(alleles_indices)
                            if '/' in alleles_indices:
                                allele1 = alleles_indices.split('/')[0]
                                allele2 = alleles_indices.split('/')[1]
                            elif '|' in alleles_indices:
                                allele1 = alleles_indices.split('|')[0]
                                allele2 = alleles_indices.split('|')[1]
                            try:
                                base1 = alleles[int(allele1)]
                            except ValueError:
                                base1 = '.'
                            try:
                                base2 = alleles[int(allele2)]
                            except ValueError:
                                base2 = '.'
                        seqs[i] += str(base1)
                        seqs[i] += str(base2)

        with open(out_file, 'w') as f:
            for i in range(len(sample_name)):
                fasta_str = f">{sample_name[i]}\n{seqs[i]}\n"
                f.writelines(fasta_str)

    def sv_tree(self):
        self.convert_vcf_to_seq(self.vcf_file, f'{self.out_path}/SV.matrix')
        infile = f'{self.out_path}/SV.matrix'
        if self.tree_software == "iqtree":
            cmd.built_tree(infile, self.out_path, "SV", self.thread)
        else:
            print(f'tree_software: {self.tree_software} can not be recognized, please check it')
            sys.exit()
