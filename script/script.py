import os
import sys
import subprocess
import multiprocessing
import re
import time

from configparser import ConfigParser
from collections import OrderedDict
from Bio import SeqIO

import script.get_parser as get_parser

opt_cfg = get_parser.get_parser()


class RunCmd():
    def __init__(self):
        self.software_path, self.opt = self.get_config()
        self.aln_software = 'mafft'
        self.tree_software = 'raxmal'
        self.thread = 10
        self.cds_mrege = None
        self.pep_mrege = None
        self.codon_mrege = None
        self.out_path = os.getcwd()
        for k, v in self.software_path.items():
            setattr(self, str(k), v)
        for k, v in self.opt.items():
            setattr(self, str(k), v)
            # print(str(k), v)
        if int(self.thread) <= 1:
            self.thread = 2
        # self.check_software() #运行时先检查软件是否可用

    def get_config(self):
        cfg_parser = ConfigParser()
        # read options
        cfg_parser.read(opt_cfg)
        opt = dict(cfg_parser.items('opt'))
        software_path = dict(cfg_parser.items('software'))

        return software_path, opt

    def check_software(self):
        """遍历配置文件中的每个软件，检查其是否可用"""
        for software_path in self.software_path.values():
            software_path = software_path.strip()  # 去除路径前后的空格
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
        """根据目录获取目录中所有文件，返回一个列表"""
        file_list = [f for f in os.listdir(in_path) if f != '']
        return [os.path.join(in_path, infile) for infile in file_list]

    def mkdir(self):
        """创建文件夹"""
        out_dir = ['01_cds_format', '02_pep', '03_hmm_out', '04_OG', '05_seq/01_cds',  '05_seq/02_pep', '06_aln/01_aln',
                   '06_aln/02_codon', '06_aln/03_trim', '06_aln/03_trim_rename', '07_tree/01_coatree',
                   '07_tree/02_contree', '08_result']
        for d in out_dir:
            try:
                os.makedirs(f"{self.out_path}/{d}")
            except OSError:
                pass

    def merge_files(self, directory):
        """将所有cds文件合并为一个文件，便于提取，返回一个字符串"""
        combined = ''
        for filename in os.listdir(directory):
            filepath = os.path.join(directory, filename)
            if os.path.isfile(filepath):
                with open(filepath, 'r') as f:
                    combined += f.read()
        return combined

    def get_seq_by_id(self, idfile, seqfile, result_file):
        """根据基因id列表提取序列，seqfile为一个字符串"""
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
        """多进程运行format_and_trans函数"""
        in_path = self.in_path
        in_file_list = self.get_infile_list(in_path)
        p = multiprocessing.Pool(int(self.thread))
        p.map(self.format_and_trans, in_file_list)
        p.close()
        p.join()

    def rename_id(self, infile, outfile):
        """将基因id仅保留物种名部分"""
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
        """运行hmmscan"""
        outfile = f"{self.out_path}/03_hmm_out/{os.path.splitext(os.path.split(infile)[1])[0]}.tbl"
        run = f'{self.hmmscan} --tblout {outfile} --noali -E 1e-50 --cpu 1 {self.orthodb} {infile} >/dev/null 2>&1'
        # run = f'echo {infile} > {outfile}'
        status = self.run_command(run)
        # print(status)
        return status

    def run_hmmscan_pl(self):
        """多进程运行hmmscan"""
        infile_list = self.get_infile_list(self.out_path + '/02_pep')
        # print(infile_list)
        p = multiprocessing.Pool(int(self.thread))
        statuss = p.map(self.HMMscan, infile_list)
        p.close()
        p.join()
        return statuss

    def aln(self, infile, outfile):
        """序列比对"""
        if self.aln_software == 'mafft':
            run = '{} {} > {} 2>/dev/null'.format(self.mafft, infile, outfile)
        else:
            run = '{} -align {} -output {}'.format(self.muscle, infile, outfile)
        status = self.run_command(run)
        return status

    def trim(self, infile, outfile):
        """序列修剪"""
        run = '{} -in {} -out {} -automated1'.format(self.trimal, infile, outfile)
        status = self.run_command(run)
        return status

    def built_tree(self, infile, outpath, basename, thread):
        """构建树"""
        if self.tree_software.upper() == 'IQTREE':
            run = f'{self.iqtree} -s {infile} -pre {outpath}/{basename} -nt {thread} -m MFP --quiet -B 1000'
        elif self.tree_software.upper() == 'FASTTREE':
            run = f'{self.fasttree} {infile} > {outpath}/{basename}.tre'
        else:
            if self.seq.upper() in ['CDS', 'CODON']:
                run = f'{self.raxmal} -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -s {basename} -n {thread} -w {outpath}'
            else:
                run = f'{self.raxmal} -f a -x 12345 -p 12345 -# 100 -m PROTGAMMAJTT -s {basename} -n {thread} -w {outpath}'
        return self.run_command(run)

    def get_genetree(self, infile):
        """串联并联法构建基因树流程，这里再加一个判断，建树序列为cds、pep或codon"""
        if self.seq.upper() == 'CDS':
            if self.cds_mrege is None:
                cds_dir = f"{self.out_path}/01_cds_format/"
                self.cds_mrege = self.merge_files(cds_dir)
        elif self.seq.upper() == 'PEP':
            if self.pep_mrege is None:
                pep_dir = f"{self.out_path}/02_pep/"
                self.pep_mrege = self.merge_files(pep_dir)
        else:
            if self.codon_mrege is None:
                codon_dir = f"{self.out_path}/06_aln/02_codon"
                self.codon_mrege = self.merge_files(codon_dir)
        basename = os.path.splitext(os.path.split(infile)[1])[0]
        seq = f"{self.out_path}/05_seq/{basename}.fa"
        aln = f"{self.out_path}/06_aln/01_aln/{basename}.aln"
        trim = f"{self.out_path}/06_aln/03_trim/{basename}.trim"
        trim_rename = f"{self.out_path}/06_aln/03_trim_rename/{basename}.trim"
        tree_path = f"{self.out_path}/07_tree/01_coatree"
        # 提取cds、比对、修剪、修改id、构建树
        self.get_seq_by_id(infile, self.cds_mrege, seq)
        self.aln(seq, aln)
        self.trim(aln, trim)
        self.rename_id(trim, trim_rename)
        self.built_tree(trim_rename, tree_path, basename, 2)

    def run_genetree_mul(self):
        """多进程构建并联基因树"""
        OG_list = self.get_infile_list(f"{self.out_path}/04_OG/")
        # print(infile_list)
        n = int(self.thread)//2
        p = multiprocessing.Pool(n)
        statuss = p.map(self.get_genetree, OG_list)
        p.close()
        p.join()
        return statuss

    def merge_gene_trees(self, inpath, outpath, soft):
        """合并基因树，这里需要注意有的文件（iqtree）直接合并会没有换行符，后面需要实际跑一下看结果"""
        out_file = f"{outpath}/merge_gene_trees.tre"
        with open(out_file, 'w') as out:
            if str(soft).upper() == "FASTTREE":
                for f in os.listdir(inpath):
                    f= f'{inpath}/{f}'
                    if os.path.splitext(f)[1] == '.tre':
                        with open(f, 'r') as infile:
                            for line in infile.readlines():
                                print(line + '\n', file=out)
            elif str(soft).upper() == "RAXMAL":
                for f in os.listdir(inpath):
                    if os.path.splitext(os.path.split(f)[1])[0] == "RAxML_bipartitions":
                        with open(f, 'r') as infile:
                            for line in infile.readlines():
                                print(line + '\n', file=out)
            else:
                for f in os.listdir(inpath):
                    if os.path.splitext(f)[1] == '.treefile':
                        with open(f, 'r') as infile:
                            for line in infile.readlines():
                                print(line + '\n', file=out)

    def run_astral(self):
        self.merge_gene_trees(f"{self.out_path}/07_tree/01_coatree", f"{self.out_path}/08_result", self.tree_software)
        run = f'java -jar {self.astral} -i {self.out_path}/08_result/merge_gene_trees.tre' \
              f' -r 100 -o {self.out_path}/08_result/Astral_tree.nwk'
        return self.run_command(run)

    def get_super_gene(self):
        """将trim_rename后的每个OG串联为super gene， 为确保每个物种序列长度一致，当某个OG中物种覆盖率不是100%时，
        用’-‘*seq_length作为缺失物种的序列。输入参数为物种基因list目录，以及序列比对后的OG目录（文件为fasta格式），
        将物种基因list目录下文件名作为串联后序列名称"""
        OG_list = self.get_infile_list(f'{self.out_path}/06_aln/02_trim_rename')
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
        """构建并联树"""
        self.get_super_gene()
        self.built_tree(f'{self.out_path}/06_aln/03_con.trim',
                        f'{self.out_path}/07_tree/02_contree', 'contree', self.thread)
