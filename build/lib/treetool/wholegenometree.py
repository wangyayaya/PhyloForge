"""whole genoeme sequence to tree"""

import os
import warnings
import shutil

from Bio import SeqIO
from Bio import BiopythonWarning
from treetool import script
from treetool.script import RunCmd
cmd = script.RunCmd()

warnings.simplefilter('ignore', BiopythonWarning)


class WholeGenomeTree(RunCmd):
    def __init__(self):
        super().__init__()

    def mkdir_(self):
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

    def aln_whole_genome(self, infile, outfile):
        """Sequence alignment"""
        """除了mafft另外两个软件没有设置线程数选项，muscle有个快捷模式，clustalw什么都没有？"""
        if str(self.aln_software) == 'muscle':
            run = '{} -align {} -output {} -maxiters 1 -diags'.format(self.muscle, infile, outfile)
        elif str(self.aln_software) == 'mafft':
            run = '{} --thread {} {} > {} 2>/dev/null'.format(self.mafft, self.thread, infile, outfile)
        elif str(self.aln_software) == 'clustalw':
            run = f'{self.clustalw} -INFILE={infile} -OUTFILE={outfile} -OUTPUT=FASTA'
        status = cmd.run_command(run)
        return status

    def run_whole_genome_tree(self):
        self.mkdir_()
        self.formatid(self.in_path, f'{self.out_path}/01_sequence_format/')
        self.aln_whole_genome(f'{self.out_path}/01_sequence_format/seq.fa', f'{self.out_path}/02_aln/01_aln/seq.aln')
        cmd.trim(f'{self.out_path}/02_aln/01_aln/seq.aln', f'{self.out_path}/02_aln/02_trim/seq.trim')
        cmd.built_tree(f'{self.out_path}/02_aln/02_trim/seq.trim', f'{self.out_path}/03_tree', 'WholeGenome', self.thread)
        seq_type = 'WholeGenome'
        if self.tree_software == "fasttree":
            shutil.copyfile(f'{self.out_path}/03_tree/{seq_type}.tre',
                            f'{self.out_path}/04_result/{seq_type}_fasttree.nwk')
        elif self.tree_software.upper() == "iqtree":
            shutil.copyfile(f'{self.out_path}/03_tree/{seq_type}.treefile',
                            f'{self.out_path}/04_result/{seq_type}_iqtree.nwk')
        else:
            shutil.copyfile(f'{self.out_path}/03_tree/RAxML_bipartitions.{seq_type}',
                            f'{self.out_path}/04_result/{seq_type}_RAxML_bipartitions.nwk')
            shutil.copyfile(f'{self.out_path}/03_tree/RAxML_bipartitionsBranchLabels.{seq_type}',
                            f'{self.out_path}/04_result/{seq_type}_RAxML_bipartitionsBranchLabels.nwk')
