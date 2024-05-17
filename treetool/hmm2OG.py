import os
import multiprocessing
import shutil
import datetime

from treetool.script import RunCmd


# 筛选hmmscan结果到OG，新加orthofinder结果到OG
class HMM_OG(RunCmd):
    def __init__(self):
        super().__init__()

        self.raw_path = f'{self.in_path}'
        self.hmm_in_path = f'{self.out_path}/03_hmm_out'
        self.og_out_path = f'{self.out_path}/04_OG'

        # 防止每次运行生成的OGs扰乱本次结果，每次都将原OGs删掉再筛选OGs
        try:
            shutil.rmtree(f'{self.og_out_path}')
        except FileNotFoundError:
            pass
        try:
            os.remove(f'{self.out_path}/map.txt')
        except FileNotFoundError:
            pass
        os.mkdir(f'{self.og_out_path}')

        try:
            self.cover = int(self.cover)
        except AttributeError:
            self.cover = len(self.get_file_paths())

        try:
            self.copy_number = int(self.copy_number)
        except AttributeError:
            self.copy_number = 1

        current_time = datetime.datetime.now()
        print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] Selecting OGs...")

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
            file_paths.append(os.path.join(self.hmm_in_path, file))
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
            out_file = f'{self.og_out_path}/busco{busco_id}'
            with open(out_file, 'w') as f:
                for i in tmp:
                    print(i, file=f)
                    # print(i)

    def run_hmm2OG(self):
        """multiprocess run save_best_match_to_file"""
        all_busco_id = self.get_all_busco_id()
        best_match_list = self.get_all_best_match()
        with multiprocessing.Pool(int(self.thread)) as pool:
            pool.starmap(self.save_best_match_to_file, [(busco_id, best_match_list) for busco_id in all_busco_id])
        if not os.listdir(self.og_out_path):
            print(f"The number of selected OGs is 0. Please adjust the cover parameter or copy_number parameter, "
                  f"change the mode parameter to 2, and rerun the task until an appropriate number of OGs is selected."
                  f"The selected OGs can be viewed in {self.og_out_path}. ")
        else:
            file_count = len(os.listdir(self.og_out_path))
            current_time = datetime.datetime.now()
            print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] A total of {file_count} OGs were filtered out. ")

    def run_hmm2OG_multi_copy(self):
        # 把所有结果转换为一个列表result：[{OG1:sp1:[gene1, gene2],sp2:[gene1, gene2]},OG2:{……},OG3:{……}]
        # 每个OG保留物种所有拷贝，可以使用Astral-pro（https://github.com/chaoszhang/A-pro）合并基因树，但是不能做串联树
        # 这一步筛选出的OG与单拷贝模式有所区别，单拷贝模式是保留busco与gene互相最佳匹配，
        # 单拷贝模式中某一物种拷贝数大于设置的参数是即将改物种删掉，这个模式中这种情况则是将整个OG删掉，所以晒出的OG可能会少于上过模式
        all_busco_list = self.get_all_busco_id()
        result = []
        for busco_id in all_busco_list:
            OG = {}
            OG[busco_id] = {}
            for in_file in self.get_file_paths():
                sp = os.path.splitext(os.path.split(in_file)[1])[0]
                OG[busco_id][sp] = []
                with open(in_file, 'r') as f:
                    for line in f:
                        if not line.startswith('#'):
                            line = line.split()
                            busco = line[0]
                            gene_id = line[2]
                            if busco == busco_id:
                                OG[busco_id][sp].append(gene_id)
                # 根据copy_number进行筛选
                if len(OG[busco_id][sp]) > int(self.copy_number) or len(OG[busco_id][sp]) == 0:
                    OG[busco_id].pop(sp)
            result.append(OG)
        # print(result)
        # 根据cover进行筛选OG至04_OG目录并保留result中符合条件的结果
        # 这里可以多进程，但是单进程跑也就十几分钟懒得改了
        for og in result:
            for k, v in og.items():
                cover = len(og[k])
                if cover >= int(self.cover):
                    # print(k)
                    # print(v.keys())
                    result.remove(og)
                    out_OG = f'{self.og_out_path}/busco{k}'
                    map_out = f'{self.out_path}/map.txt'
                    with open(out_OG, 'w') as f1, open(map_out, 'a') as f2:
                        for gene in v.values():
                            for gene in gene:
                                species_name = gene.split('|')[0][0:]
                                print(gene, file=f1)  # 保留合格的至04_OG目录
                                # 保留基因及其对应物种名，后续astral-pro合并需要用到
                                ### astral-pro直接合并基因树问题太多了……
                                print(f'{gene}\t{species_name}', file=f2)
        if not os.listdir(self.og_out_path):
            print(f"The number of selected OGs is 0. Please adjust the cover parameter or copy_number parameter, "
                  f"change the mode parameter to 2, and rerun the task until an appropriate number of OGs is selected."
                  f"The selected OGs can be viewed in {self.og_out_path}.")
        else:
            file_count = len(os.listdir(self.og_out_path))
            current_time = datetime.datetime.now()
            print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] A total of {file_count} OGs were filtered out.")
        return result

    def ortho2OG(self, retain_mult_copy):
        or_path = f'{self.out_path}/03_orthofinder/OrthoFinder'
        for item in os.listdir(or_path):
            if os.path.isdir(os.path.join(or_path, item)) and item.startswith('Results_'):
                or_result_dir = f'{or_path}/{item}'
                # print(or_result_dir)

        with open(f"{or_result_dir}/Orthogroups/Orthogroups.tsv") as f:
            next(f)
            for l in f.readlines():
                cover_tmp = 0
                line = l.split('\t')
                for gene in line[1:]:
                    gene = gene.strip().split(',')
                    if len(gene) <= self.copy_number:
                        if len(gene) == 1:
                            if gene[0] == '':
                                pass
                            else:
                                cover_tmp += 1
                        else:
                            cover_tmp += 1

                # 保留多拷贝基因做树这样就算完成了，保留单拷贝的话就保留第一个
                if cover_tmp >= self.cover:
                    with open(f'{self.og_out_path}/{line[0]}', 'w') as fw:
                        for gene in line[1:]:
                            gene = gene.strip().split(',')
                            if len(gene) <= self.copy_number:
                                if gene[0] == '':
                                    pass
                                else:
                                    if retain_mult_copy:
                                        for id in gene:
                                            print(id.strip(), file=fw)
                                    else:
                                        print(gene[0], file=fw)
        if not os.listdir(self.og_out_path):
            print(f"The number of selected OGs is 0. Please adjust the cover parameter or copy_number parameter, "
                  f"change the mode parameter to 2, and rerun the task until an appropriate number of OGs is selected."
                  f"The selected OGs can be viewed in {self.og_out_path}.")
        else:
            file_count = len(os.listdir(self.og_out_path))
            current_time = datetime.datetime.now()
            print(f"[{current_time.strftime('%Y-%m-%d %H:%M:%S')}] A total of {file_count} OGs were filtered out.")
