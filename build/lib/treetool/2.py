"""snp vcf to tree"""
from itertools import islice
def convert_vcf_to_seq(vcf_file):
    with open(vcf_file) as vcf_file:
        s = ''
        n = 0
        while True:
            next_n_lines = list(islice(vcf_file, 163540))
            if not next_n_lines:
                break
            # 处理读取到的行数据
            for line in next_n_lines:
                if line.startswith("#CHROM"):
                    sample_name = line.strip().split()[9:]
                    print(sample_name)
                    seqs = [""] * len(sample_name)
                    # break
                    # print(seqs)
                else:
                    s += line[0]
                    n += 1
                    print(line[0], n)
    print(s)

convert_vcf_to_seq("..\茶树泛基因组\TBF.vcf")
print('2,islice')