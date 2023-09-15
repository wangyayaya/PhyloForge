import os
import multiprocessing
import sys
import shutil

from configparser import ConfigParser


from est import run



class HMM_OG:
    def __init__(self):
        self.thread = 10
        opt = self.get_config()
        for k, v in opt.items():
            setattr(self, str(k), v)

        self.raw_path = f'{self.in_path}'
        self.in_path = f'{self.out_path}/03_hmm_out'
        self.out_path = f'{self.out_path}/04_OG'

        # 防止每次运行生成的OGs扰乱本次结果，每次都将原OGs删掉再筛选OGs
        shutil.rmtree(f'{self.out_path}')
        os.mkdir(f'{self.out_path}')

        try:
            cover = int(self.cover)
        except AttributeError:
            self.cover = len(self.get_file_paths())

        try:
            copy_number = int(self.copy_number)
        except AttributeError:
            self.copy_number = 1

    def get_config(self):
        opt_cfg = run.get_parser()[0]
        cfg_parser = ConfigParser()
        # read options
        cfg_parser.read(opt_cfg)
        opt = dict(cfg_parser.items('lcn_opt'))
        return opt

    def get_file_paths(self):
        """这里输入文件改成与原输入文件对应，这样就可以在原输入目录中增删文件实现续跑"""
        file_list = [f for f in os.listdir(self.raw_path)]
        base_name_list = []
        for f in file_list:
            base_name = os.path.splitext(os.path.split(f)[1])[0]
            base_name_list.append(base_name)
        # print(base_name_list)

        file_paths = []
        for file in base_name_list:
            file = f'{file}.tbl'
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
            # 选择单拷贝或低拷贝，1：单拷贝，0以上低拷贝
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
        """save the result by get_best_match as list"""
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
        """multiprocess run save_best_match_to_file"""
        all_busco_id = self.get_all_busco_id()
        best_match_list = self.get_all_best_match()
        with multiprocessing.Pool(int(self.thread)) as pool:
            pool.starmap(self.save_best_match_to_file, [(busco_id, best_match_list) for busco_id in all_busco_id])
