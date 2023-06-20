import os
import sys
import subprocess
import multiprocessing
import re
import time
import shutil

from configparser import ConfigParser
from collections import OrderedDict
from Bio import SeqIO

import get_parser as get_parser

opt_cfg = get_parser.get_parser()


class RunCmd():
    def __init__(self):
        self.aln_software = 'mafft'
        self.tree_software = 'raxmal'
        self.seq = 'cds'
        self.thread = 10
        self.cds_mrege = None
        self.pep_mrege = None
        self.codon_mrege = None
        self.out_path = os.getcwd()

        self.software_path, self.opt = self.get_config()
        for k, v in self.software_path.items():
            setattr(self, str(k), v)
        for k, v in self.opt.items():
            setattr(self, str(k), v)
            # print(str(k), v)

        if int(self.thread) <= 1:
            self.thread = 2

        self.check_software()  # 运行时先检查软件是否可用

    def get_config(self):
        cfg_parser = ConfigParser()
        # read options
        cfg_parser.read(opt_cfg)
        opt = dict(cfg_parser.items('opt'))
        software_path = dict(cfg_parser.items('software'))

        return software_path, opt

    def check_software(self):
        """Iterate through each software in the profile to check if it is available"""
        for software_path in self.software_path.values():
            software_path = software_path.strip()
            if not os.path.exists(software_path):
                try:
                    subprocess.run([software_path, "--version"], check=True,
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # 执行软件，检查执行结果
                except FileNotFoundError:
                    print(f"Error: {software_path} is not available!")
                    exit(1)
                except Exception:
                    continue

    def get_infile_list(self, in_path):
        """Gets all the files in the directory according to the directory and returns a list"""
        file_list = [f for f in os.listdir(in_path) if f != '']
        return [os.path.join(in_path, infile) for infile in file_list]

    def mkdir(self):
        """make directory"""
        out_dir = ['01_cds_format', '02_pep', '03_hmm_out', '04_OG', '05_seq/01_cds', '05_seq/02_pep', '06_aln/01_aln',
                   '06_aln/02_trim', '06_aln/03_trim_rename', '07_tree/01_coatree', '07_tree/02_contree', '08_result']
        for d in out_dir:
            try:
                os.makedirs(f"{self.out_path}/{d}")
            except OSError:
                pass

    def merge_files(self, directory):
        """Merge all files in the directory into one file for easy extraction and return a string"""
        combined = ''
        for filename in os.listdir(directory):
            filepath = os.path.join(directory, filename)
            if os.path.isfile(filepath):
                with open(filepath, 'r') as f:
                    combined += f.read()
        return combined

    def get_seq_by_id(self, idfile, seqfile, result_file):
        """The sequence is extracted from the list of gene ids, seqfile as a string format"""
        aDict = OrderedDict()
        for line in seqfile.splitlines():
            if line.startswith('>'):
                k = line.split()[0][1:]  # 去掉开头的'>',只取序列的名字
                aDict[k] = []
            else:
                aDict[k].append(line)
        # print(aDict)

        with open(idfile) as idf:
            with open(result_file, 'w') as result:
                for line in idf:
                    name = line.strip()
                    if aDict[name]:
                        print(">%s\n%s" % (name, ''.join(aDict[name])), file=result)

    def format_and_trans(self, in_file):
        """将用户提供原始cds序列格式化基因id即翻译为蛋白序列，新id为物种名（basename infile）_cds_raw_id,
        后续以_cds_为分隔符保留物种名作为基因id构建基因树，便于合并基因树为物种树"""
        base_name = os.path.splitext(os.path.split(in_file)[1])[0]
        cds_fm_out = f"{self.out_path}/01_cds_format/{base_name}.cds"
        cds_fm_out = open(cds_fm_out, 'w')
        pep_out = f"{self.out_path}/02_pep/{base_name}.pep"
        pep_out = open(pep_out, 'w')
        for line in SeqIO.parse(in_file, 'fasta'):
            raw_id = line.id
            new_id = f'>{base_name}_cds_{raw_id}'
            cds = re.sub(r'[^ATCGUatcgu]', 'N', str(line.seq))  # 用正则表达式将非法字符替换为N
            pep = line.seq.translate(table="Standard")
            print(f'{new_id} \n {cds}', file=cds_fm_out)
            print(f'{new_id} \n {pep}', file=pep_out)

    # 多进程
    def run_format_and_trans(self):
        """Run the format_and_trans function in multiple processes"""
        in_path = self.in_path
        in_file_list = self.get_infile_list(in_path)
        p = multiprocessing.Pool(int(self.thread))
        p.map(self.format_and_trans, in_file_list)
        p.close()
        p.join()

    def rename_id(self, infile, outfile):
        """Keep only the species name part of the gene id in the sequence"""
        with open(infile, 'r') as f1, open(outfile, 'w') as f2:
            for line in f1:
                if line.startswith('>'):
                    species_name = line.split('_cds_')[0][1:]
                    f2.write('>' + species_name + '\n')
                else:
                    f2.write(line)

    def run_command(self, command):
        try:
            output = subprocess.check_output(command, shell=True, stderr=subprocess.DEVNULL)
            return 0, output.decode()
        except subprocess.CalledProcessError as e:
            print(e)
            return e.returncode

    def HMMscan(self, infile):
        """Run hmmscan"""
        outfile = f"{self.out_path}/03_hmm_out/{os.path.splitext(os.path.split(infile)[1])[0]}.tbl"
        run = f'{self.hmmscan} --tblout {outfile} --noali -E 1e-50 --cpu 1 {self.orthodb} {infile} >/dev/null 2>&1'
        # run = f'echo {infile} > {outfile}'
        status = self.run_command(run)
        # print(status)
        return status

    def run_hmmscan_pl(self):
        """Run hmmscan in multiple processes"""
        self.mkdir()
        self.run_format_and_trans()
        infile_list = self.get_infile_list(self.out_path + '/02_pep')
        # print(infile_list)
        p = multiprocessing.Pool(int(self.thread))
        statuss = p.map(self.HMMscan, infile_list)
        p.close()
        p.join()
        return statuss

    def aln(self, infile, outfile):
        """Sequence alignment"""
        if self.aln_software == 'mafft':
            run = '{} {} > {} 2>/dev/null'.format(self.mafft, infile, outfile)
        else:
            run = '{} -align {} -output {}'.format(self.muscle, infile, outfile)
        status = self.run_command(run)
        return status

    def aln_codon(self, alnfile, cdsfile, out):
        """Convert protein sequence alignment to codon"""
        run = f"perl {self.pal2nal} {alnfile} {cdsfile} -output fasta >{out}"
        return self.run_command(run)

    def trim(self, infile, outfile):
        """Trim sequence"""
        run = '{} -in {} -out {} -automated1'.format(self.trimal, infile, outfile)
        status = self.run_command(run)
        return status

    def built_tree(self, infile, outpath, basename, thread):
        """Construct tree"""
        if self.tree_software.upper() == 'IQTREE':
            run = f'{self.iqtree} -s {infile} -pre {outpath}/{basename} -nt {thread} -m MFP --quiet -B 1000'
        elif self.tree_software.upper() == 'FASTTREE':
            run = f'{self.fasttree} {infile} > {outpath}/{basename}.tre'
        else:
            if self.seq.upper() in ['CDS', 'CODON']:
                run = f'{self.raxml} -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -s {infile} -n {basename} -T {thread} -w {outpath}'
            else:
                run = f'{self.raxml} -f a -x 12345 -p 12345 -# 100 -m PROTGAMMAJTT -s {infile} -n {basename} -T {thread} -w {outpath}'
        return self.run_command(run)

    def get_genetree(self, infile):
        """pipeline to Construct gene tree"""
        def tree_cds_pep(seq_type):
            basename = os.path.splitext(os.path.split(infile)[1])[0]
            basename = f"{basename}_{seq_type}"
            seq = f"{self.out_path}/05_seq/{basename}.fa"
            aln = f"{self.out_path}/06_aln/01_aln/{basename}.aln"
            trim = f"{self.out_path}/06_aln/02_trim/{basename}.trim"
            trim_rename = f"{self.out_path}/06_aln/03_trim_rename/{basename}.trim"
            tree_path = f"{self.out_path}/07_tree/01_coatree"
            # 提取cds、比对、修剪、修改id、构建树
            if seq_type == 'cds':
                self.get_seq_by_id(infile, self.cds_mrege, seq)
            elif seq_type == 'pep':
                self.get_seq_by_id(infile, self.pep_mrege, seq)
            self.aln(seq, aln)
            self.trim(aln, trim)
            self.rename_id(trim, trim_rename)
            self.built_tree(trim_rename, tree_path, basename, 2)

        def tree_codon():
            basename = os.path.splitext(os.path.split(infile)[1])[0]
            seq_cds = f"{self.out_path}/05_seq/{basename}_cds.fa"
            seq_pep = f"{self.out_path}/05_seq/{basename}_pep.fa"
            aln_pep = f"{self.out_path}/06_aln/01_aln/{basename}_pep.aln"
            codon_seq = f"{self.out_path}/06_aln/01_aln/{basename}_codon.aln"
            trim = f"{self.out_path}/06_aln/02_trim/{basename}_codon.trim"
            trim_rename = f"{self.out_path}/06_aln/02_trim_rename/{basename}_codon.trim"
            tree_path = f"{self.out_path}/07_tree/01_coatree"
            # 提取cds、比对、修剪、修改id、构建树
            self.get_seq_by_id(infile, self.cds_mrege, seq_cds)
            self.get_seq_by_id(infile, self.pep_mrege, seq_pep)
            self.aln(seq_pep, aln_pep)
            self.aln_codon(aln_pep, seq_cds, codon_seq)
            self.trim(codon_seq, trim)
            self.rename_id(trim, trim_rename)
            self.built_tree(trim_rename, tree_path, f"{basename}_codon", 2)

        if self.seq.upper() == 'CDS':
            if self.cds_mrege is None:
                cds_dir = f"{self.out_path}/01_cds_format/"
                self.cds_mrege = self.merge_files(cds_dir)
            tree_cds_pep("cds")
        elif self.seq.upper() == 'PEP':
            if self.pep_mrege is None:
                pep_dir = f"{self.out_path}/02_pep/"
                self.pep_mrege = self.merge_files(pep_dir)
            tree_cds_pep('pep')
        else:
            if self.cds_mrege or self.pep_mrege is None:
                pep_dir = f"{self.out_path}/02_pep/"
                self.pep_mrege = self.merge_files(pep_dir)
                cds_dir = f"{self.out_path}/01_cds_format/"
                self.cds_mrege = self.merge_files(cds_dir)
            tree_codon()

    def run_genetree_mul(self):
        """Multiprocess construct trees"""
        OG_list = self.get_infile_list(f"{self.out_path}/04_OG/")
        # print(infile_list)
        n = int(self.thread) // 2
        p = multiprocessing.Pool(n)
        statuss = p.map(self.get_genetree, OG_list)
        p.close()
        p.join()
        return statuss

    def merge_gene_trees(self, inpath, outpath, soft):
        """Merge gen trees"""
        out_file = f"{outpath}/merge_gene_trees.tre"
        with open(out_file, 'w') as out:
            if str(soft).upper() == "FASTTREE":
                for f in os.listdir(inpath):
                    f = f'{inpath}/{f}'
                    if os.path.splitext(f)[1] == '.tre':
                        with open(f, 'r') as infile:
                            for line in infile.readlines():
                                print(line + '\n', file=out)
            elif str(soft).upper() == "RAXML":
                for f in os.listdir(inpath):
                    f = f'{inpath}/{f}'
                    if os.path.splitext(os.path.split(f)[1])[0] == "RAxML_bipartitions":
                        with open(f, 'r') as infile:
                            for line in infile.readlines():
                                print(line + '\n', file=out)
            else:
                for f in os.listdir(inpath):
                    f = f'{inpath}/{f}'
                    if os.path.splitext(f)[1] == '.treefile':
                        with open(f, 'r') as infile:
                            for line in infile.readlines():
                                print(line + '\n', file=out)

    def run_astral(self):
        """Construct Coalescence tree"""
        self.merge_gene_trees(f"{self.out_path}/07_tree/01_coatree", f"{self.out_path}/08_result", self.tree_software)
        run = f'java -jar {self.astral} -i {self.out_path}/08_result/merge_gene_trees.tre' \
              f' -r 100 -o {self.out_path}/08_result/Astral_tree.nwk'
        return self.run_command(run)

    def get_super_gene(self):
        """将trim_rename后的每个OG串联为super gene， 为确保每个物种序列长度一致，当某个OG中物种覆盖率不是100%时，
        用’-‘*seq_length作为缺失物种的序列。将物种基因list目录下文件名作为串联后序列名称"""
        OG_list = self.get_infile_list(f'{self.out_path}/06_aln/03_trim_rename')
        sp_list = self.get_infile_list(self.in_path)
        result = open(f'{self.out_path}/06_aln/03_con.trim', 'w')

        for sp_path in sp_list:
            sp = os.path.splitext(os.path.split(sp_path)[1])[0]
            gene_list = [sp]

            sp_seq = ''
            for OG in OG_list:
                seq_tmp = ''
                # id_in_genelist = False
                for line in SeqIO.parse(OG, 'fasta'):
                    gene_id = line.id
                    # print(id)
                    seq = line.seq
                    seq_len = len(seq)
                    # while not id_in_genelist:
                    if gene_id in gene_list:
                        seq_tmp = seq
                        # id_in_genelist = True
                        break
                    else:
                        seq_tmp = '-' * seq_len
                        # break
                sp_seq += seq_tmp
            print(f'>{sp}\n{sp_seq}', file=result)

    def run_contree(self):
        """Construct Concatenation tree"""
        self.get_super_gene()
        self.built_tree(f'{self.out_path}/06_aln/03_con.trim',
                        f'{self.out_path}/07_tree/02_contree', 'contree', self.thread)
        if self.tree_software.upper() == "FASTTREE":
            shutil.copyfile(f'{self.out_path}/07_tree/02_contree/contree.tre', f'{self.out_path}/08_result/contree_fasttree.nwk')

        elif self.tree_software.upper() == "RAXML":
            shutil.copyfile(f'{self.out_path}/07_tree/02_contree/RAxML_bipartitions.contree',
                            f'{self.out_path}/08_result/contree_RAxML_bipartitions.nwk')
            shutil.copyfile(f'{self.out_path}/07_tree/02_contree/RAxML_bipartitionsBranchLabels.contree',
                            f'{self.out_path}/08_result/contree_RAxML_bipartitionsBranchLabels.nwk')
        else:
            shutil.copyfile(f'{self.out_path}/07_tree/02_contree/contree.treefile',
                            f'{self.out_path}/08_result/contree_iqtree.nwk')
