import os
import sys
import re
import warnings
import shutil

from Bio import SeqIO
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

from treetool import script, run
cmd = script.RunCmd()


class OrganelleTree():
    def __init__(self):
        self.coa_con = 0
        opt_cfg = run.get_parser()[0]
        opt = cmd.get_config(opt_cfg, 'organelle_opt')
        for k, v in opt.items():
            setattr(self, str(k), v)

        self.check()

    def check(self):
        try:
            int(self.cover)
        except AttributeError:
            pass
        except ValueError:
            print('Please check the cover, which must be a int')
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
            if str(self.seq).upper() not in ["CDS", "PEP", "CODON"]:
                print(f'The seq: {self.seq} is not recognized, and CDS will be used')
        except AttributeError:
            print('The seq is not specified, and CDS will be used')

        try:
            int(self.thread)
        except AttributeError:
            pass
        except ValueError:
            print("thread should be int (>=2), please check it")
            sys.exit()

        try:
            if str(self.aln_software).upper() not in ['MAFFT', 'MUSCLE']:
                print(f'The aln_software: {self.aln_software} is not recognized, and MAFFT will be used or sequence alignment')
        except AttributeError:
            print('Sequence alignment software is not specified, MAFFT will be used for sequence alignment')

        try:
            if str(self.tree_software).upper() not in ['RAXML', 'IQTREE', 'FASTTREE']:
                print(f'The tree_software: {self.tree_software} is not recognized, and RAxML will be used')
        except AttributeError:
            print('The tree_software is not specified, raxml software will be used')


    def mkdir(self):
        """make directory"""
        out_dir = ['01_cds_format', '02_pep', '03_OG_all', '04_OG', '05_seq/', '06_aln/01_aln', '06_aln/02_trim',
                   '06_aln/03_trim_rename', '07_tree/01_coatree', '07_tree/02_contree', '08_result']
        for d in out_dir:
            try:
                # shutil.rmtree(f'{self.out_path}/{d}')
                os.makedirs(f'{self.out_path}/{d}')
            except OSError:
                pass

    def trans_format_getgenelist(self, inpath, cds_out_path, pep_out_path, OG_path):
        fa_dict = {}
        in_file_list = [f for f in os.listdir(inpath) if f != '']
        in_file_path = [os.path.join(inpath, infile) for infile in in_file_list]
        sp_list = [f.split('.')[0] for f in in_file_list]
        gene_list = []
        # trans and format
        for f in in_file_path:
            sp_gene_list = []
            base_name = os.path.splitext(os.path.split(f)[1])[0]
            cds_out = open(f'{cds_out_path}/{base_name}.cds', 'w')
            pep_out = open(f'{pep_out_path}/{base_name}.pep', 'w')
            for line in SeqIO.parse(f, 'fasta'):
                try:
                    if line.description.split('[')[1].split(']')[0].split('=')[0].upper() == 'GENE':
                        gene_id = line.description.split('[')[1].split(']')[0].split('=')[1].upper()
                    else:
                        gene_id = line.id
                except IndexError:
                    gene_id = line.id
                OG = open(f'{OG_path}/{gene_id}', 'a')
                new_id = f'{base_name}_wy_{gene_id}'
                if gene_id not in sp_gene_list:
                    sp_gene_list.append(gene_id)
                    print(new_id, file=OG)
                if gene_id not in gene_list:
                    gene_list.append(gene_id)
                cds = re.sub(r'[^ATCGUatcgu]', 'N', str(line.seq))  # 这里感觉也没必要
                pep = line.seq.translate(table="Standard")
                print(f'>{new_id} \n {cds}', file=cds_out)
                print(f'>{new_id} \n {pep}', file=pep_out)
            # print(sp_gene_list)
        # print(gene_list)
        return gene_list

    def select_OG(self):
        for OG in self.trans_format_getgenelist(self.in_path, f'{self.out_path}/01_cds_format/',
                                                f'{self.out_path}/02_pep', f'{self.out_path}/03_OG_all'):
            count = len(open(f'{self.out_path}/03_OG_all/{OG}', 'r').readlines())
            if count >= int(self.cover):
                shutil.copyfile(f'{self.out_path}/03_OG_all/{OG}', f'{self.out_path}/04_OG/{OG}')

    def run_organelle_tree(self):
        self.mkdir()
        self.select_OG()
        cmd.run_genetree_mul()
        cmd.run_astral()
        cmd.run_contree()
