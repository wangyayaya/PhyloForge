import os
import sys
import multiprocessing

from configparser import ConfigParser

import get_parser as get_parser
import script as sccript

opt_cfg = get_parser.get_parser()


class HMM_OG:
    def __init__(self):
        self.thread = 10
        opt = self.get_config()
        for k, v in opt.items():
            setattr(self, str(k), v)

        self.in_path = f'{self.out_path}/03_hmm_out'
        self.out_path = f'{self.out_path}/04_OG'

        try:
            self.cover
        except AttributeError:
            self.cover = len(self.get_file_paths())

    def get_config(self):
        cfg_parser = ConfigParser()
        # read options
        cfg_parser.read(opt_cfg)
        opt = dict(cfg_parser.items('opt'))
        return opt

    def get_file_paths(self):
        """根据输入目录获取所有文件绝对路径，保存为一个列表，输入目录作为一个参数或直接读取self.in_path？"""
        file_paths = []
        for file in os.listdir(self.in_path):
            file_paths.append(os.path.join(self.in_path, file))
        return file_paths

    def get_all_busco_id(self):
        """将所有物种hmm结果的busco id提取为一个列表"""
        all_busco_list = []
        for file_path in self.get_file_paths():
            with open(file_path, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        busco_id = line.split()[0]
                        if busco_id not in all_busco_list:
                            all_busco_list.append(busco_id)
        # print(all_busco_list)
        return all_busco_list

    def get_best_match(self, in_file):
        """获取每个文件的基因与busco id分别获取最佳匹配，类似双向blast，返回字典，key：busco id；value：gene id"""
        gene_busco_dict = {}
        busco_gene_dict = {}
        busco_list = {}
        gene_list = {}
        with open(in_file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    line = line.split()
                    busco_id = line[0]
                    gene_id = line[2]
                    e_true = float(eval(line[4]))
                    if busco_id not in busco_list:
                        busco_list[busco_id] = {}
                    busco_list[busco_id][gene_id] = e_true
                    if gene_id not in gene_list:
                        gene_list[gene_id] = {}
                    gene_list[gene_id][busco_id] = e_true
        # print(gene_list)
        # print(busco_list)
        for gene_id, busco_dict in gene_list.items():
            e_tmp = 1
            for busco_id, e_true in busco_dict.items():
                if e_true < e_tmp:
                    gene_busco_dict[gene_id] = busco_id
                    e_tmp = e_true
        for busco_id, gene_dict in busco_list.items():
            # 选择单拷贝或低拷贝，1：单拷贝，0-3或0-2……低拷贝
            if len(gene_dict) <= int(self.copy_number):
                e_tmp = 1
                for gene_id, e_true in gene_dict.items():
                    if e_true < e_tmp:
                        busco_gene_dict[busco_id] = gene_id
                        e_tmp = e_true
        # with open(f'{self.tmp_path}/{os.path.splitext(os.path.split(in_file)[1])[0]}.tmp', 'w') as f:
        best_match = {}
        for k, v in gene_busco_dict.items():
            if v in busco_gene_dict and busco_gene_dict[v] == k:
                best_match[v] = k
                # print(v)
                # print(f'{v}\t{k}', file=f)
        # print(best_match.keys())
        return best_match

    def get_all_best_match(self):
        """将get_best_match函数得到的所有物种最佳匹配的结果以列表形式储存"""
        best_match_list = []
        for in_file in self.get_file_paths():
            busco_gene_best = self.get_best_match(in_file)
            best_match_list.append(busco_gene_best)
        return best_match_list

    def save_best_match_to_file(self, busco_id, best_match_list):
        tmp = []
        for best_match in best_match_list:
            for k, v in best_match.items():
                # print(k)
                if k == busco_id:
                    tmp.append(v)
        if len(tmp) >= int(self.cover):
            out_file = f'{self.out_path}/busco{busco_id}'
            with open(out_file, 'w') as f:
                for i in tmp:
                    print(i, file=f)
                    #print(i)

    def run_hmm2OG(self):
        """多进程运行save_best_match_to_file函数"""
        all_busco_id = self.get_all_busco_id()
        best_match_list = self.get_all_best_match()
        with multiprocessing.Pool(int(self.thread)) as pool:
            pool.starmap(self.save_best_match_to_file, [(busco_id, best_match_list) for busco_id in all_busco_id])
