"""snp vcf to tree"""

def convert_vcf_to_seq(vcf_file):
    with open(vcf_file) as vcf_file:
        s = ''
        n = 0
        for line in vcf_file:
            if line.startswith("#CHROM"):
                sample_name = line.strip().split()[9:]
                print(sample_name)
                seqs = [""] * len(sample_name)
                #break
                # print(seqs)
            else:
                s += line[0]
                n += 1
                print(line[0], n)
    print(s)

convert_vcf_to_seq("..\茶树泛基因组\TBF.vcf")
print('test,正常')