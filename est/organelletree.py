import os
import re
import warnings
import shutil

from Bio import SeqIO
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

from est import script, run

cmd = script.RunCmd()


class OrganelleTree():
    def __init__(self):
        self.coa_con = 0
        self.mode = 0
        opt_cfg = run.get_parser()[0]
        opt = cmd.get_config(opt_cfg, 'organelle_opt')
        for k, v in opt.items():
            setattr(self, str(k), v)

    def mkdir(self):
        """make directory"""
        out_dir = ['01_cds_format', '02_pep', '03_OG_all', '04_OG', '05_seq/', '06_aln/01_aln', '06_aln/02_trim',
                   '06_aln/03_trim_rename', '07_tree/01_coatree', '07_tree/02_contree', '08_result']
        for d in out_dir:
            try:
                shutil.rmtree(f'{self.out_path}/{d}')
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
                    gene_id = line.description.split('[')[1].split(']')[0].split('=')[1].upper()
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
