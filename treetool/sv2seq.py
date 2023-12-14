#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Function: Convert SNP data in VCF format to sequence (fasta/phylip)
# wang 20231108

import argparse
import sys


def convert_vcf_to_seq(vcf_file, out_file):
    """只对插入和缺失进行转换，插入时：
    如果存在杂合突变：
    ref=0, alt=1, [0,1]
    0/0: 00
    0/1: 01
    1/1: 11
    如果没有杂合突变：
    0/0: 0
    1/1: 1
    缺失时：
    ref=1, alt=0, [1,0]
    如果存在杂合突变：
    0/0: 11
    0/1: 10
    1/1: 00
    如果没有杂合突变：
    0/0: 1
    1/1: 0
    """

    with open(vcf_file) as vcf_file:
        lines = vcf_file.readlines()
    for line in lines:
        if line.startswith("#CHROM"):
            sample_name = line.strip().split()[9:]
            # print(sample_name)
            seqs = [""] * len(sample_name)
            # print(seqs)
        elif line.startswith("#"):
            pass
        elif line.strip():
            fields = line.strip().split()
            sv_type = fields[4]
            if sv_type == '<DEL>':
                alleles = [1, 0]
            elif sv_type == '<INS>':
                alleles = [0, 1]
            # print(alleles)
            # 这里好恶心
            # 虽然很恶心，但是改了又报错或者转换的结果不对
            # 恶心就恶心吧，不改了
            Hom = False
            for i, alleles_indices in enumerate(fields[9:], start=0):
                if alleles_indices.split(":")[0]:
                    alleles_indices = alleles_indices.split(":")[0]
                    # print(alleles_indices)
                    try:
                        if '/' in alleles_indices:
                            allele1 = alleles_indices.split('/')[0]
                            allele2 = alleles_indices.split('/')[1]
                        elif '|' in alleles_indices:
                            allele1 = alleles_indices.split('|')[0]
                            allele2 = alleles_indices.split('|')[1]
                        base1 = alleles[int(allele1)]
                        base2 = alleles[int(allele2)]
                    except ValueError:
                        base1 = base2 = '.'
                    if allele1 == allele2:
                        Hom = True

            if Hom:
                for i, alleles_indices in enumerate(fields[9:], start=0):
                    if alleles_indices.split(":")[0]:
                        alleles_indices = alleles_indices.split(":")[0]
                        # print(alleles_indices)
                        try:
                            if '/' in alleles_indices:
                                allele1 = alleles_indices.split('/')[0]
                                allele2 = alleles_indices.split('/')[1]
                            elif '|' in alleles_indices:
                                allele1 = alleles_indices.split('|')[0]
                                allele2 = alleles_indices.split('|')[1]
                            base = alleles[int(allele1)]
                        except ValueError:
                            base = '.'
                    seqs[i] += str(base)
            else:
                for i, alleles_indices in enumerate(fields[9:], start=0):
                    if alleles_indices.split(":")[0]:
                        alleles_indices = alleles_indices.split(":")[0]
                        # print(alleles_indices)
                        if '/' in alleles_indices:
                            allele1 = alleles_indices.split('/')[0]
                            allele2 = alleles_indices.split('/')[1]
                        elif '|' in alleles_indices:
                            allele1 = alleles_indices.split('|')[0]
                            allele2 = alleles_indices.split('|')[1]
                        base1 = alleles[int(allele1)]
                        base2 = alleles[int(allele2)]
                    seqs[i] += str(base1)
                    seqs[i] += str(base2)

    with open(out_file, 'w') as f:
        for i in range(len(sample_name)):
            fasta_str = f">{sample_name[i]}\n{seqs[i]}\n"
            f.writelines(fasta_str)







convert_vcf_to_seq("AMP_SV_head1000.vcf", "tmp1.fa")