# EasySpeciesTree使用说明

## 1.所需文件

该软件需要两个需要用户提供的文件：

1）需要做物种树的物种cds序列（建议每个基因只保留最长转录本），存放于一个文件夹下，物种cds文件建议以genus_species.cds (fa/fas/fasta)来命名，最终物种树中即以genus_species作为物种名。

2）与做物种树的物种对应的BUSCO单拷贝数据库，该文件可以按如下步骤获取：

1，访问BUSCO数据库：[Index of /v5/data/lineages/ (ezlab.org)](https://busco-data.ezlab.org/v5/data/lineages/)或者https://busco-archive.ezlab.org/v3/ 下载与物种对应的数据库，例如被子植物一般下载有胚植物（Embryophyta odb10）单拷贝基因数据集：

```shell
$ wget -c https://busco-archive.ezlab.org/v3/datasets/prerelease/embryophyta_odb10.tar.gz -O embryophyta_odb10.tar.gz
$ tar -zxvf embryophyta_odb10.tar.gz
$ cd embryophyta_odb10
$ cat hmms/*.hmm >embryophyta_odb10.hmm
$ hmmpress embryophyta_odb10.hmm
# 获取embryophyta_odb10.hmm路径，将其写入对应的配置文件中
$ realpath embryophyta_odb10.hmm
```

## 2.所需软件

1）hmmsearch软件，安装该软件后，主要使用hmmpress命令和hmmscan命令

2）序列比对软件，可选择mafft或muscle，默认为mafft

3）序列修剪软件：trimal

4）构树软件：可选择fasttree、iqtree或raxmal，默认为raxmal

5）基因树合并软件：astral

## 3.设置配置文件

可通过run.py -g >run.config获取配置文件

run.config参数说明:

```
[opt]
in_path = D:\WANG\Species_Tree\cs\
out_path = D:\WANG\Species_Tree\Species_Tree\tmp\
thread = int(>=2)
orthodb = /path/to/buscoDB
aln_software = muscle/mafft（）
tree_software = iqtree/fasttree/raxmal(Multithreaded version)
cover = 7 #每个OG中物种数量，当做树物种为100，该参数选择50时，即每个OG中物种覆盖率至少为50%
copy_number = 1/2/3…… #基因拷贝数，1：单拷贝，2以上为低拷贝
step = 0 #选择运行步骤，
#0：从头开始跑完全流程
#从头跑到筛选完OG（初步判断OG数目是否适当）
#单独运行筛选OG这一步骤（筛选出适当的OG数目）
#从筛选完OG后运行完后续全步骤
seq = cds/pep/codon #选择用那种类型序列来构建物种树
coa_con = 0/1/2 #选择构建并联树或串联树

[software]
hmmscan=/path/to/hmmscan
mafft=/path/to/mafft
muscle=/path/to/muscle
trimal=/path/to/trimal
raxmal=/path/to/raxmal
iqtree=/path/to/iqtree
astral=/path/to/astral.jar
```

运行程序：run.py -c run.config

#### 4.结果文件说明

