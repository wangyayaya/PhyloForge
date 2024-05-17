"""snp vcf to tree"""
import re
import os
import shutil
import datetime
import sys

from Bio import SeqIO
import multiprocessing
from treetool import script
from treetool.script import RunCmd
cmd = script.RunCmd()


class VcfTree(RunCmd):
    def __init__(self):
        super().__init__()
        try:
            os.makedirs(f"{self.out_path}/working_dir")
            self.working_dir = os.path.abspath(f"{self.out_path}/working_dir")
        except OSError:
            shutil.rmtree(f"{self.out_path}/working_dir")
            os.makedirs(f"{self.out_path}/working_dir")
            self.working_dir = os.path.abspath(f"{self.out_path}/working_dir")

    def generate_ambiguous_code(self, base1, base2):
        # 非纯和突变碱基的歧义码IUPAC Ambiguity Codes，类型没有vcf2phylip的多
        # 但运行结果检查了一下和那个脚本是一样的，因为只考虑单碱基突变的话就只有下面这几种突变
        # 一个自己的想法，如果一个位点参考是A，突变是T，那么存在杂合突变时是否可以将0/0转换为AA，1/1转换为TT，而0/1转换为AT，感觉可行.
        # SNP感觉几乎每一个位点都会有杂合，所以这样无形中相当于将数据量增大了一倍！！！
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

    def cutFile(self, vcf_file, num, out_path):
        """VCF文件过大时，将文件拆分多个小文件"""
        with open(vcf_file) as f:
            for line in f:
                if line.startswith("#CHROM"):
                    header = line
                    break

        sourceFileData = open(vcf_file, 'r', encoding='utf-8')
        ListOfLine = [line for line in sourceFileData.read().splitlines() if not line.startswith("#")]
        n = len(ListOfLine)
        p = n // num + 1
        split_files = []
        for i in range(num):
            split_file = os.path.join(out_path, os.path.basename(os.path.splitext(vcf_file)[0]) + str(i) +
                                      os.path.splitext(vcf_file)[-1])
            split_files.append(split_file)
            # print(split_file)
            destFileData = open(split_file, "w", encoding='utf-8')
            destFileData.write(header)
            if (i == num - 1):
                for line in ListOfLine[i * p:]:
                    destFileData.write(line + '\n')
            else:
                for line in ListOfLine[i * p:(i + 1) * p]:
                    destFileData.write(line + '\n')
            destFileData.close()
        return split_files

    def convert_vcf_to_seq(self, vcf_file):
        with open(vcf_file) as f:
            for line in f:
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
                                    allele1, allele2 = map(int, re.split('/|\|', alleles_indices))
                                    base1 = alleles[allele1]
                                    base2 = alleles[allele2]
                                except ValueError:
                                    base1 = base2 = '.'
                                # print(base1, base2)
                                base = self.generate_ambiguous_code(base1, base2)
                                # print(base)
                                seqs[i] = "".join([seqs[i], base])

        fasta_str = ""
        for i in range(len(sample_name)):
            fasta_str += f">{sample_name[i]}\n{seqs[i]}\n"
        out_str = fasta_str
        out_file = os.path.splitext(vcf_file)[0] + '.fa'
        #print(out_file)
        with open(out_file, 'w') as f:
            f.writelines(out_str)
        return out_file

    def merge_seq(self, input_files, out_file, format):
        """merge split fa files"""
        sequences = {}
        for file_name in input_files:
            for seq_record in SeqIO.parse(file_name, "fasta"):
                seq_id = seq_record.id
                # 判断当前 ID 是否已经存在于 sequences 字典中
                if seq_id in sequences:
                    # 如果存在，则将当前序列拼接到已有序列的末尾
                    sequences[seq_id] += seq_record.seq
                else:
                    # 如果不存在，则将当前序列添加到字典中
                    sequences[seq_id] = seq_record.seq

        if format == 'phylip':
            out_f = out_file + '.phy'
            with open(out_f, 'w') as file:
                num_sequences = len(sequences)
                sequence_length = len(next(iter(sequences.values())))
                file.write(f'{num_sequences} {sequence_length}\n')
                for sequence_id, sequence in sequences.items():
                    file.write(f'{sequence_id}  {sequence}\n')
        else:
            out_f = out_file + '.fa'
            with open(out_f, "w") as f:
                for seq_id, seq in sequences.items():
                    f.write(">{}\n{}\n".format(seq_id, seq))



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

    def change_file_extension(self, file_list, new_extension):
        new_file_list = []
        for file_name in file_list:
            base_name, _ = os.path.splitext(file_name)
            new_file_name = base_name + new_extension
            new_file_list.append(new_file_name)
        return new_file_list

    def snp_tree(self):
        current_time = datetime.datetime.now()
        print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] Converting a VCF file to a base sequence...")
        basename = 'SNP'
        splitfiles = self.cutFile(self.vcf_file, int(self.thread), self.working_dir)
        fa_list = self.change_file_extension(splitfiles, '.fa')
        try:
            p = multiprocessing.Pool(int(self.thread))
            p.map(self.convert_vcf_to_seq, splitfiles)
            current_time = datetime.datetime.now()
            print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] Finish Converting a VCF file to a base sequence.")
        except Exception:
            current_time = datetime.datetime.now()
            print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] Failed to Convert VCF file to base sequence, "
                  f"please check the input file!")
            sys.exit(1)

        if self.tree_software.upper() == 'TREEBEST' or self.tree_software.upper() == 'FASTTREE':
            self.merge_seq(fa_list, f'{self.out_path}/all_snp', 'fa')
            infile = f'{self.out_path}/all_snp.fa'
        else:
            self.merge_seq(fa_list, f'{self.out_path}/all_snp', 'phylip')
            infile = f'{self.out_path}/all_snp.phy'
            print(infile)

        current_time = datetime.datetime.now()
        print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] Constructing tree...")
        status = cmd.built_tree(infile, f'{self.out_path}', basename, 1)
        if status == 0:
            print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] Finish constructing the SNP tree.")
        else:
            print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] SNP tree failed to construct.")



