import os
import subprocess
import multiprocessing
import re
import shutil
import sys
import datetime
import glob
import warnings
from configparser import ConfigParser
from collections import OrderedDict
from Bio import SeqIO
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

from treetool import run


class RunCmd:
    def __init__(self):
        self.aln_software = 'mafft'
        self.tree_software = 'raxml'
        self.seq = 'cds'
        self.thread = 10
        self.cds_mrege = None
        self.pep_mrege = None
        self.codon_mrege = None
        self.out_path = os.getcwd()

        opt_cfg, tree = run.get_parser()[0], run.get_parser()[1]
        if not os.path.isfile(opt_cfg):
            print(f"{opt_cfg} is not a valid file, please check it")
            print("Exiting...")
            sys.exit(1)

        self.software_path = self.get_soft_path()
        for k, v in self.software_path.items():
            setattr(self, str(k).lower(), v)

        # 把所有基于低拷贝的合并在一起
        if tree == 'lcn':
            self.opt = self.get_config(opt_cfg, 'lcn_opt')
        elif tree == 'organelle':
            self.opt = self.get_config(opt_cfg, 'organelle_opt')
        elif tree == 'snp':
            self.opt = self.get_config(opt_cfg, 'snp_opt')
        elif tree == 'whole_genome':
            self.opt = self.get_config(opt_cfg, 'whole_genome_opt')
        elif tree == 'sv':
            self.opt = self.get_config(opt_cfg, 'sv_opt')
            self.tree_software = 'iqtree'
        elif tree == 'gene':
            self.opt = self.get_config(opt_cfg, 'gene_opt')
            self.retain_multi_copy = True

        for k, v in self.opt.items():
            if 'software' in k:
                setattr(self, str(k), v.lower())
            else:
                setattr(self, str(k), v)

        try:
            os.makedirs(f"{self.out_path}")
            self.out_path = os.path.abspath(f"{self.out_path}")
        except OSError:
            self.out_path = os.path.abspath(f"{self.out_path}")

        try:
            if int(self.thread) <= 1:
                self.thread = 2
        except ValueError:
            print("thread should be int, please check it")
            sys.exit()
        try:
            if str(self.tree_software) not in ['raxml', 'iqtree', 'fasttree', 'treebest', 'phyml']:
                print(f'The tree_software: {self.tree_software} is not recognized, please check it')
                sys.exit()
        except AttributeError:
            print('The tree_software is not specified, raxml software will be used for lcn/organelle/whole_genome'
                  'iqtree software will be used for sv, treebest software will be used for snp')

        if tree == 'lcn':
            try:
                if str(self.retain_multi_copy).upper() == 'TRUE':
                    self.retain_multi_copy = True
                    self.coa_con = 0
                    try:
                        os.path.abspath(os.path.dirname(self.astral_pro))
                    except AttributeError:
                        current_time = datetime.datetime.now()
                        print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] If [retain_multi_copy = Ture], astral-pro "
                              f"software must be configured, please refer to the README file for configuration.")
                        sys.exit(1)
                else:
                    self.retain_multi_copy = False
            except AttributeError:
                self.retain_multi_copy = False

    def get_config(self, opt_cfg, opt):
        cfg_parser = ConfigParser()
        # read options
        cfg_parser.read(opt_cfg)
        opt = dict(cfg_parser.items(opt))
        return opt

    def get_soft_path(self):
        cfg_parser = ConfigParser()
        ob_path = os.path.split(__file__)[0]
        cfg = os.path.join(ob_path, "cfg", "software.ini")
        cfg_parser.read(cfg)
        software_path = dict(cfg_parser.items('software'))
        return software_path

    def check_software(self):
        """Iterate through each software in the profile to check if it is available"""
        """这里用conda打包了，可以不用这个函数了"""
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
        out_dir = ['01_cds_format', '02_pep', '03_hmm_out', '04_OG', '05_seq/', '06_aln/01_aln',
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
                        print(">%s\n%s" % (name, ''.join(aDict[name]).rstrip()), file=result)

    def format_and_trans(self, in_file):
        """The original CDS sequences provided by the user will be formatted, with gene IDs translated into protein
        sequences. The new IDs will be in the format of species name (basename infile)|raw_id. Subsequently,
        the gene tree will be constructed using the species name preserved as the gene ID, using | as a separator.
        This facilitates the merging of gene trees into a species tree."""
        base_name = os.path.splitext(os.path.split(in_file)[1])[0]
        cds_fm_out = f"{self.out_path}/01_cds_format/{base_name}.cds"
        cds_fm_out = open(cds_fm_out, 'w')
        pep_out = f"{self.out_path}/02_pep/{base_name}.pep"
        pep_out = open(pep_out, 'w')
        for line in SeqIO.parse(in_file, 'fasta'):
            try:
                raw_id = line.id
                new_id = f'>{base_name}|{raw_id}'.replace(':', '_')  # 基因ID中不能有冒号
                cds = re.sub(r'[^ATCGUatcgu]', 'N', str(line.seq))  # 将非法字符替换为N
                pep = line.seq.translate(table="Standard")
                # 重新添加orthofinder后不能识别有些异常字符，所以这里需要再对pep再格式化 *暂未
                print(f'{new_id}\n{cds}', file=cds_fm_out)
                print(f'{new_id}\n{pep}', file=pep_out)
            except Exception:
                print(f"Please check the input file: {in_file}, the sequence may have abnormal characters ")
                pass

    # 多进程
    def run_format_and_trans(self):
        """Run the format_and_trans function in multiple processes"""
        in_path = self.in_path
        in_file_list = self.get_infile_list(in_path)
        file_count = len(os.listdir(in_path))
        current_time = datetime.datetime.now()
        print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] Formatting and translating {file_count} CDS files...")
        try:
            p = multiprocessing.Pool(int(self.thread))
            p.map(self.format_and_trans, in_file_list)
            p.close()
            p.join()
            current_time = datetime.datetime.now()
            print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] Finish formatting and translating all CDS files.")
        except Exception:
            print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] Failed to format and translate CDS files.")

    def rename_id(self, infile, outfile):
        """Keep only the species name part of the gene id in the sequence"""
        with open(infile, 'r') as f1, open(outfile, 'w') as f2:
            for line in f1:
                if line.startswith('>'):
                    species_name = line.split('|')[0][1:]
                    f2.write('>' + species_name + '\n')
                else:
                    f2.write(line)

    def run_command(self, command):
        """上传conda后，明明hmmscan运行是成功的，还是会返回非0，而且也不报错，非常奇怪"""
        try:
            subprocess.check_output(command, shell=True, stderr=subprocess.DEVNULL)
            return 0
        except subprocess.CalledProcessError as e:
            print(e)
            sys.exit(e.returncode)

    def HMMscan(self, infile):
        """Run hmmscan"""
        outfile = f"{self.out_path}/03_hmm_out/{os.path.splitext(os.path.split(infile)[1])[0]}.tbl"
        run = f'{self.hmmscan} --tblout {outfile} --noali -E 1e-50 --cpu 1 {self.orthodb} {infile}'

        if os.path.exists(outfile):
            with open(outfile, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    line = line.strip()
                    if "[ok]" in line:
                        status = 0
                        break
                else:
                    status = self.run_command(run)
        else:
            # run = f'echo {infile} > {outfile}'
            status = self.run_command(run)
        if status != 0:
            print(f'The file {infile} failed to run hmmsearch program!')
        return status

    def run_hmmscan_pl(self):
        """Run hmmscan in multiple processes"""
        self.mkdir()
        self.run_format_and_trans()
        current_time = datetime.datetime.now()
        print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] Runing hmmscan program, this may take a long time...")
        infile_list = self.get_infile_list(self.out_path + '/02_pep')
        # print(infile_list)
        p = multiprocessing.Pool(int(self.thread))
        statuss = p.map(self.HMMscan, infile_list)
        p.close()
        p.join()
        end_time = datetime.datetime.now()
        if any(status != 0 for status in statuss):
            print(f"[{end_time.strftime('%Y-%m-%d %H:%M:%S')}] The hmmscan program is finished."
                  f"The above file fails to run. Please check and restart.")
            sys.exit(1)
        else:
            print(f"[{end_time.strftime('%Y-%m-%d %H:%M:%S')}] The hmmscan program is successfully finished.")
        return statuss

    def run_orthofinder(self):
        # 只要跑orthofinder就将原文件全删掉，防止生成多个orthofinder的结果目录
        self.mkdir()
        self.run_format_and_trans()
        current_time = datetime.datetime.now()
        print(
            f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] Runing OrthoFinder program, this may take a long time...")
        try:
            shutil.rmtree(f'{self.out_path}/03_orthofinder')
        except Exception:
            pass
        shutil.copytree(f"{self.out_path}/02_pep", f"{self.out_path}/03_orthofinder")
        orthofinder_dir = f"{self.out_path}/03_orthofinder"
        orthofinder_cmd = f'orthofinder -t {int(self.thread)} -f {orthofinder_dir} -S diamond -og'
        status = self.run_command(orthofinder_cmd)
        end_time = datetime.datetime.now()
        if status != 0:
            print(f"{end_time.strftime('%Y-%m-%d %H:%M:%S')} The OrthoFinder program failed to run.")
            sys.exit(1)
        else:
            print(f"[{end_time.strftime('%Y-%m-%d %H:%M:%S')}] The OrthoFinder program is successfully finished.")
        return status

    def aln(self, infile, outfile):
        """Sequence alignment"""
        if str(self.aln_software) == 'muscle':
            run = '{} -align {} -output {}'.format(self.muscle, infile, outfile)
        elif str(self.aln_software) == 'mafft':
            run = '{} {} > {}'.format(self.mafft, infile, outfile)
        elif str(self.aln_software) == 'clustalw':
            run = f'{self.clustalw} -INFILE={infile} -OUTFILE={outfile} -OUTPUT=FASTA'
        status = self.run_command(run)
        return status

    def aln_codon(self, alnfile, cdsfile, out):
        """Convert protein sequence alignment to codon"""
        run = f"{self.pal2nal} {alnfile} {cdsfile} -output fasta >{out}"
        # 提取密码子第几位 codon:不移除，codon1保留1，2位，condon2：保留第三位
        status = self.run_command(run)
        return status

    def codon_position_select(self, infile, outfile, pos):
        with open(infile) as f, open(outfile, 'w') as outfile:
            for line in SeqIO.parse(f, 'fasta'):
                codon = line.seq
                outfile.writelines(f'>{line.id}\n')
                seq = ''
                if pos == 'codon1':
                    for i in range(0, len(codon), 3):
                        # print(i)
                        pos1 = i
                        pos2 = i + 1
                        # print(pos2)
                        seq += codon[pos1]
                        seq += codon[pos2]
                    outfile.writelines(f'{seq}\n')
                elif pos == 'codon2':
                    for i in range(2, len(codon), 3):
                        seq += codon[i]
                    outfile.writelines(f'{seq}\n')
                else:
                    sys.exit()

    def trim(self, infile, outfile):
        """Trim sequence"""
        run = '{} -in {} -out {} -automated1'.format(self.trimal, infile, outfile)
        status = self.run_command(run)
        return status

    def built_tree(self, infile, outpath, basename, thread):
        """Construct tree"""
        if os.path.isfile(infile) and os.path.getsize(infile) == 0:
            pass
        else:
            if self.tree_software == 'iqtree':
                run = f'{self.iqtree} -s {infile} -pre {outpath}/{basename} -nt {thread} -m MFP --quiet -B 1000 -redo'
            elif self.tree_software == 'fasttree':
                run = f'{self.fasttree} {infile} > {outpath}/{basename}.tre'
            elif self.tree_software == 'raxml':
                if self.seq in ['cds', 'codon', 'codon1', 'codon2']:
                    run = f'{self.raxml} -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -s {infile} -n {basename} -T {thread} -w {outpath}'
                else:
                    run = f'{self.raxml} -f a -x 12345 -p 12345 -# 100 -m PROTGAMMAJTT -s {infile} -n {basename} -T {thread} -w {outpath}'
            elif self.tree_software == 'treebest':
                run = f'{self.treebest} nj -b 1000  {infile} >{outpath}/{basename}_treebest.NHX'
            elif self.tree_software == 'phyml':
                run = f'{self.phyml} -i {infile} -b 100 -m HKY85 -f m -v e -a e -o tlr'
            status = self.run_command(run)
        return status

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
            # extract cds sequence、alignment、trim、rename id、construct tree
            if seq_type == 'cds':
                self.get_seq_by_id(infile, self.cds_mrege, seq)
            elif seq_type == 'pep':
                self.get_seq_by_id(infile, self.pep_mrege, seq)
            self.aln(seq, aln)
            self.trim(aln, trim)
            if not self.retain_multi_copy:
                self.rename_id(trim, trim_rename)
                self.built_tree(trim_rename, tree_path, basename, 2)
            else:
                self.built_tree(trim, tree_path, basename, 2)

        def tree_codon():
            basename = os.path.splitext(os.path.split(infile)[1])[0]
            seq_cds = f"{self.out_path}/05_seq/{basename}_cds.fa"
            seq_pep = f"{self.out_path}/05_seq/{basename}_pep.fa"
            aln_pep = f"{self.out_path}/06_aln/01_aln/{basename}_pep.aln"
            codon_seq = f"{self.out_path}/06_aln/01_aln/{basename}_codon.aln"
            trim = f"{self.out_path}/06_aln/02_trim/{basename}_{self.seq}.trim"
            trim_rename = f"{self.out_path}/06_aln/03_trim_rename/{basename}_{self.seq}.trim"
            tree_path = f"{self.out_path}/07_tree/01_coatree"
            # 提取cds、比对、修剪、修改id、构建树
            self.get_seq_by_id(infile, self.cds_mrege, seq_cds)
            self.get_seq_by_id(infile, self.pep_mrege, seq_pep)
            self.aln(seq_pep, aln_pep)
            self.aln_codon(aln_pep, seq_cds, codon_seq)
            if os.path.isfile(codon_seq) and os.path.getsize(codon_seq) == 0:
                """to codon， len of cds is not mul of 3. codon file is empty, trimal failed"""
                os.remove(codon_seq)
            else:
                # 在这里筛选codon位点
                codon_seq_select = f"{self.out_path}/06_aln/01_aln/{basename}_{self.seq}.aln"
                if self.seq.upper() in ['CODON1', 'CODON2']:
                    self.codon_position_select(codon_seq, codon_seq_select, self.seq.lower())
                    codon_seq = codon_seq_select
                else:
                    pass

                self.trim(codon_seq, trim)
                if not self.retain_multi_copy:
                    self.rename_id(trim, trim_rename)
                    self.built_tree(trim_rename, tree_path, f"{basename}_{self.seq}", 2)
                else:
                    self.built_tree(trim, tree_path, f"{basename}_{self.seq}", 2)

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
        """一开始只考虑做lcn，这个脚本中很多函数不太好利用"""
        current_time = datetime.datetime.now()
        print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] Constructing trees...")

        OG_list = self.get_infile_list(f"{self.out_path}/04_OG/")
        # 删除过去跑出的树及比对后后序列，防止续跑时前面的结果影响
        dir_list = ['05_seq/', '06_aln/01_aln', '06_aln/02_trim', '06_aln/03_trim_rename',
                    '07_tree/01_coatree', '07_tree/02_contree']
        for dir in dir_list:
            try:
                shutil.rmtree(f'{self.out_path}/{dir}')
            except FileNotFoundError:
                pass
            os.mkdir(f'{self.out_path}/{dir}')

        if len(OG_list) <= 0:
            print("The number of OGs is 0, exit......")
            sys.exit()
        else:
            # print(infile_list)

            if self.tree_software.upper() == 'FASTTREE':
                n = int(self.thread)
            else:
                n = int(self.thread) // 2

            p = multiprocessing.Pool(n)
            statuss = p.map(self.get_genetree, OG_list)
            p.close()
            p.join()
            return statuss

    def merge_gene_trees(self, inpath, outpath, soft):
        """Merge gen trees"""
        out_file = f"{outpath}/merge_gene_trees_{self.seq.lower()}_{soft}.tre"

        if str(soft).upper() == "FASTTREE":
            gene_tree_list = glob.glob(f'{inpath}/*_{self.seq.lower()}.tre')
        elif str(soft).upper() == "IQTREE":
            gene_tree_list = glob.glob(f'{inpath}/*_{self.seq.lower()}.treefile')
        else:
            gene_tree_list = glob.glob(f'{inpath}/RAxML_bipartitions.*_{self.seq.lower()}')

        with open(out_file, 'w') as out:
            for f in gene_tree_list:
                with open(f, 'r') as infile:
                    for line in infile.readlines():
                        print(line + '\n', file=out)

        return out_file

    def run_astral(self):
        """Construct Coalescence tree"""
        merge_gene_trees = self.merge_gene_trees(f"{self.out_path}/07_tree/01_coatree", f"{self.out_path}/08_result",
                                                 self.tree_software)
        if not self.retain_multi_copy:
            run = f'{self.astral} -i {merge_gene_trees} -r 100 -o ' \
                  f'{self.out_path}/08_result/coalescent-based_{self.seq.lower()}_{self.tree_software}.nwk'
        else:
            absp_astral = os.path.abspath(os.path.dirname(self.astral_pro))
            astrallib = os.path.join(absp_astral, 'lib')
            run = f'java -Djava.library.path={astrallib} -jar {self.astral_pro} -i {merge_gene_trees} ' \
                  f'-o {self.out_path}/08_result/coalescent-based_{self.seq.lower()}_{self.tree_software}_mcl.nwk ' \
                  f'-a {self.out_path}/map.txt'
        status = self.run_command(run)
        if status == 0:
            current_time = datetime.datetime.now()
            print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] Finish constructing the coalescent-based tree.")
        else:
            print(status)

    def get_super_gene(self):
        """将trim_rename后的每个OG串联为super gene， 为确保每个物种序列长度一致，当某个OG中物种覆盖率不是100%时，
        用’-‘*seq_length作为缺失物种的序列。将物种基因list目录下文件名作为串联后序列名称"""
        sp_list = self.get_infile_list(self.in_path)
        OG_list = glob.glob(f'{self.out_path}/06_aln/02_trim/*_{self.seq.lower()}.trim')
        supergene = f'{self.out_path}/06_aln/supergene_{self.seq.lower()}.trim'
        seq_type = f'concatenated_{self.seq.lower()}'

        result = open(supergene, 'w')

        for sp_path in sp_list:
            sp = os.path.splitext(os.path.split(sp_path)[1])[0]
            gene_list = [sp]

            sp_seq = ''
            for OG in OG_list:
                seq_tmp = ''
                # id_in_genelist = False
                for line in SeqIO.parse(OG, 'fasta'):
                    gene_id = line.id.split('|')[0][0:]
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
            seq_check = sp_seq.replace('-', '')
            if len(seq_check) > 0:
                print(f'>{sp}\n{sp_seq}', file=result)
            # print(f'>{sp}\n{sp_seq}')
        return supergene, seq_type

    def run_contree(self):
        """Construct Concatenation tree"""
        supergene, seq_type = self.get_super_gene()
        status = self.built_tree(supergene, f'{self.out_path}/07_tree/02_contree', seq_type, self.thread)
        if self.tree_software.upper() == "FASTTREE":
            shutil.copyfile(f'{self.out_path}/07_tree/02_contree/{seq_type}.tre',
                            f'{self.out_path}/08_result/{seq_type}_fasttree.nwk')
        elif self.tree_software.upper() == "IQTREE":

            shutil.copyfile(f'{self.out_path}/07_tree/02_contree/{seq_type}.treefile',
                            f'{self.out_path}/08_result/{seq_type}_iqtree.nwk')
        else:
            shutil.copyfile(f'{self.out_path}/07_tree/02_contree/RAxML_bipartitions.{seq_type}',
                            f'{self.out_path}/08_result/{seq_type}_RAxML_bipartitions.nwk')
            shutil.copyfile(f'{self.out_path}/07_tree/02_contree/RAxML_bipartitionsBranchLabels.{seq_type}',
                            f'{self.out_path}/08_result/{seq_type}_RAxML_bipartitionsBranchLabels.nwk')
        if status == 0:
            current_time = datetime.datetime.now()
            print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] Finish constructing the concatenation tree.")
