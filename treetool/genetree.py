import os
import shutil
import re
import multiprocessing
from Bio import SeqIO
from Bio import BiopythonWarning
import warnings

warnings.simplefilter('ignore', BiopythonWarning)

from treetool import script
from treetool.script import RunCmd

cmd = script.RunCmd()


class GeneTree(RunCmd):
    def __init__(self):
        super().__init__()

    def mkdir_(self):
        out_dir = ['01_cds_format', '02_pep', '04_OG', '05_seq/', '06_aln/01_aln', '06_aln/02_trim', '07_tree/01_coatree']
        new_dir = ["01_cds", "02_pep", "03_ID", "04_seq", "05_aln", "06_tree"]
        all_dir = out_dir + new_dir
        for d in all_dir:
            try:
                shutil.rmtree(f'{self.out_path}/{d}')
            except OSError:
                pass
            except FileNotFoundError:
                pass
        for d in out_dir:
            os.makedirs(f'{self.out_path}/{d}')

    def trans(self, in_file):
        base_name = os.path.splitext(os.path.split(in_file)[1])[0]
        cds_fm_out = f"{self.out_path}/01_cds_format/{base_name}.cds"
        cds_fm_out = open(cds_fm_out, 'w')
        pep_out = f"{self.out_path}/02_pep/{base_name}.pep"
        pep_out = open(pep_out, 'w')
        og_out = f"{self.out_path}/04_OG/{base_name}"
        og_out = open(og_out, 'w')
        for line in SeqIO.parse(in_file, 'fasta'):
            raw_id = line.id
            cds = re.sub(r'[^ATCGUatcgu]', 'N', str(line.seq))
            pep = line.seq.translate(table="Standard")
            print(f'>{raw_id}\n{cds}', file=cds_fm_out)
            print(f'>{raw_id}\n{pep}', file=pep_out)
            print(f'{raw_id}', file=og_out)

    def run_trans(self):
        """Run the format_and_trans function in multiple processes"""
        in_path = self.in_path
        in_file_list = self.get_infile_list(in_path)
        p = multiprocessing.Pool(int(self.thread))
        p.map(self.trans, in_file_list)
        p.close()
        p.join()

    def run_gene_tree(self):
        self.mkdir_()
        self.run_trans()
        cmd.run_genetree_mul()
        shutil.move(f'{self.out_path}/01_cds_format', f'{self.out_path}/01_cds')
        shutil.move(f'{self.out_path}/04_OG', f'{self.out_path}/03_ID')
        shutil.move(f'{self.out_path}/05_seq', f'{self.out_path}/04_seq')
        shutil.move(f'{self.out_path}/06_aln', f'{self.out_path}/05_aln')
        shutil.move(f'{self.out_path}/07_tree/01_coatree', f'{self.out_path}/06_tree')
        shutil.rmtree(f'{self.out_path}/07_tree')

