import sys

args = sys.argv
if len(args) <= 3:
    print(f'Usage: python {args[0]} vcf_file out_file convert_type[phy/fa]')
    sys.exit()

vcf_file = args[1]
out = args[2]
f = args[3]


def convert_vcf_to_seq(vcf_file, out, format):
    # Open the input vcf file.
    with open(vcf_file) as vcf_file:
        lines = vcf_file.readlines()

    # Extract sequence data from vcf file.
    ids = []
    seqs = []
    ids_found = False
    for line in lines:
        if line.startswith("#CHROM"):
            ids = line.strip().split()[9:]
            seqs = [""] * len(ids)
            ids_found = True
        elif ids_found and line.strip():
            fields = line.strip().split()
            ref = fields[3]
            alleles = [ref] + fields[4].split(",")
            for i, gt in enumerate(fields[9:], start=0):
                if gt.split(":")[0] == "./." or gt.split(":")[0] == ".|.":
                    gt_seq = "-"
                else:
                    if '/' in gt:
                        alleles_indices = [int(idx) for idx in gt.split("/") if idx.isdigit()]
                    elif '|' in gt:
                        alleles_indices = [int(idx) for idx in gt.split("|") if idx.isdigit()]
                    gt_seq = "".join(alleles[idx] for idx in alleles_indices if idx < len(alleles))
                    for nu in gt_seq:
                        if nu not in ['A','T','G','C','a','t','g','c','-']:
                            gt_seq = gt_seq.replace(nu, '-')
                # VCF4.2格式中，同一位点突变可能会涉及多个碱基，以最长的那一个样本的碱基长度作为所有样本长度，不足的用'-'补齐
                len_gt = 1
                for allele in alleles:
                    if len(allele) > len_gt:
                        len_gt = len(allele)
                sup_lenth = len_gt - len(gt_seq)
                gt_seq = gt_seq + '-' * sup_lenth

                seqs[i] += gt_seq

    if f == 'phy':
        # Prepare the phylip string.
        phylip_str = f"{len(ids)} {len(seqs[0])}\n"
        for i in range(len(ids)):
            phylip_str += f"{ids[i].ljust(10)} {seqs[i]}\n"
        out_str = phylip_str
    # 60个字符换一行
    else:
        fasta_str = ""
        for i in range(len(ids)):
            fasta_str += f">{ids[i]}\n"
            seq = seqs[i]
            while len(seq) > 0:
                fasta_str += seq[:60] + "\n"
                seq = seq[60:]  # 每一次删掉序列前60个碱基
        out_str = fasta_str
    # 不换行
    """else:
        fasta_str = ""
        for i in range(len(ids)):
            fasta_str += f">{ids[i]}\n{seqs[i]}\n"
        out_str = fasta_str"""

    # Write the phylip string to file.
    with open(out, "w") as out:
        out.write(out_str)

    # print(out_str)


if __name__ == '__main__':
    convert_vcf_to_seq(vcf_file, out, f)
