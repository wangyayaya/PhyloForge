"""whole genoeme sequence to tree"""

import os
import sys
import warnings
import shutil

from Bio import SeqIO
from Bio import BiopythonWarning
from treetool import script, run

cmd = script.RunCmd()

warnings.simplefilter('ignore', BiopythonWarning)


class WholeGenomeTree():
    def __init__(self):
        self.coa_con = 0
        self.mode = 0
        opt_cfg = run.get_parser()[0]
        opt = cmd.get_config(opt_cfg, 'whole_genome_opt')
        for k, v in opt.items():
            setattr(self, str(k), v)

        self.software_path = cmd.get_soft_path()
        for k, v in self.software_path.items():
            setattr(self, str(k), v)

        self.out_path = os.path.abspath(self.out_path)
        self.check()

    def check(self):
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
        out_dir = ['01_sequence_format', '02_aln/01_aln', '02_aln/02_trim', '03_tree', '04_result']
        for d in out_dir:
            try:
                os.makedirs(f'{self.out_path}/{d}')
            except OSError:
                pass

    def formatid(self, inpath, outpath):
        in_file_list = [f for f in os.listdir(inpath) if f != '']
        in_file_path = [os.path.join(inpath, infile) for infile in in_file_list]
        outfile = open(f'{outpath}/seq.fa', 'w')
        for f in in_file_path:
            base_name = os.path.splitext(os.path.split(f)[1])[0]
            for line in SeqIO.parse(f, 'fasta'):
                print(f'>{base_name}', file=outfile)
                print(line.seq, file=outfile)

    def aln(self, infile, outfile):
        """Sequence alignment"""
        if str(self.aln_software).upper() == 'MUSCLE':
            run = '{} -align {} -output {} -maxiters 1 -diags'.format(self.muscle, infile, outfile)
        else:
            run = '{} --thread {} {} > {} 2>/dev/null'.format(self.mafft, self.thread, infile, outfile)
        status = cmd.run_command(run)
        return status

    def run_whole_genome_tree(self):
        self.mkdir()
        self.formatid(self.in_path, f'{self.out_path}/01_sequence_format/')
        self.aln(f'{self.out_path}/01_sequence_format/seq.fa', f'{self.out_path}/02_aln/01_aln/seq.aln')
        cmd.trim(f'{self.out_path}/02_aln/01_aln/seq.aln', f'{self.out_path}/02_aln/02_trim/seq.trim')
        cmd.built_tree(f'{self.out_path}/02_aln/02_trim/seq.trim', f'{self.out_path}/03_tree', 'WholeGenome', self.thread)
        seq_type = 'WholeGenome'
        if self.tree_software.upper() == "FASTTREE":
            shutil.copyfile(f'{self.out_path}/03_tree/{seq_type}.tre',
                            f'{self.out_path}/04_result/{seq_type}_fasttree.nwk')
        elif self.tree_software.upper() == "IQTREE":
            shutil.copyfile(f'{self.out_path}/03_tree/{seq_type}.treefile',
                            f'{self.out_path}/04_result/{seq_type}_iqtree.nwk')
        else:
            shutil.copyfile(f'{self.out_path}/03_tree/RAxML_bipartitions.{seq_type}',
                            f'{self.out_path}/04_result/{seq_type}_RAxML_bipartitions.nwk')
            shutil.copyfile(f'{self.out_path}/03_tree/RAxML_bipartitionsBranchLabels.{seq_type}',
                            f'{self.out_path}/04_result/{seq_type}_RAxML_bipartitionsBranchLabels.nwk')
