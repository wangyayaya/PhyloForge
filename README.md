# EasySpeciesTree使用说明

## 1.所需文件

该软件需要两个需要用户提供的文件：

1）需要做物种树的物种cds序列（建议每个基因只保留最长转录本），存放于一个文件夹下，物种cds文件建议以genus_species.cds (fa/fas/fasta)来命名，最终物种树中即以genus_species作为物种名。

2）与做物种树的物种对应的BUSCO单拷贝数据库，该文件可以按如下步骤获取：

访问BUSCO数据库：[Index of /v5/data/lineages/ (ezlab.org)](https://busco-data.ezlab.org/v5/data/lineages/)或者https://busco-archive.ezlab.org/v3/ 下载与物种对应的数据库，例如被子植物一般下载有胚植物（Embryophyta odb10）单拷贝基因数据集：

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

1）hmmer软件，安装该软件后，主要使用hmmpress命令和hmmscan命令。目前大部分window版本的hmmer软件不能识别我们需要的hmm文件。**必需**

2）序列比对软件，可选择mafft或muscle，默认为mafft

3）序列修剪软件：trimal，**必需**

4）构树软件：可选择fasttree、iqtree或raxmal，默认为raxmal，最好是多线程版本，否则线程数设置不起作用

5）pal2nal，将蛋白比对序列转换为codon后做树，在配置文件中seq=codon时，需配置该软件

6）基因树合并软件：astral，**必需**

## 3.设置配置文件

### 首先配置相关软件，可运行run.py -gsc >software.config生成软件配置文件

software.config：

```
[software]
hmmscan=hmmscan
mafft=mafft
muscle=muscle
trimal=trimal
raxml=raxmlHPC-PTHREADS
iqtree=iqtree
fasttree=fasttree
astral=/opt/service/ASTRAL/astral.5.7.3.jar
pal2nal=/opt/service/miniconda3/envs/r4_py37_env/bin/pal2nal.pl
```

**注意事项：**
**#软件在环境变量中直接写调用该软件时的软件名即可，否则需要绝对/相对路径**
**#opt选项需要用到的软件必需要有对应软件，如opt中序列比对软件选择mafft，则必需配置mafft软件路径，muscle就无需配置**
**#如果选择使用codon序列构树，则软件需配置pal2nal软件**

运行run.py  -sc software.config命令配置软件

### 运行参数配置

可通过run.py -grc >run.config获取配置文件

run.config参数说明:

```
[opt]
in_path = D:\WANG\Species_Tree\cs\ #cds文件目录
out_path = D:\WANG\Species_Tree\Species_Tree\out\  #结果输出目录
thread = int(>=2) #线程数，默认10，最小2
orthodb = /path/to/buscoDB  #busco数据库路径，参考1.2
aln_software = muscle/mafft #序列比对软件，默认mafft
tree_software = iqtree/fasttree/raxmal(Multithreaded version) #构树软件，默认raxmal 
cover = 7 #每个OG中物种数量，当做树物种为100，该参数选择50时，即每个OG中物种覆盖率至少为50%。默认为物种数量，即物种覆盖率100%
copy_number = 1/2/3…… #基因拷贝数，1：单拷贝，2以上为低拷贝，默认1
mode = 0 #选择运行模式，默认0
#0：从头开始跑完全流程
#1：从头跑到筛选完OG（初步判断OG数目是否适当）
#2：单独运行筛选OG这一步骤（筛选出适当的OG数目）
#3：从筛选完OG后运行完后续全步骤
seq = cds/pep/codon #选择用那种类型序列来构建物种树，默认codon
#注意，当选择使用codon做树时，并不是所有序列都能转换为codon，因此实际用来做树的OG数目可能会少于筛选到的数目，可使用ls out_path/06_aln/02_trim|wc -l查看实际做树的OG数目
coa_con = 0/1/2 #选择构建并联树或串联树，0：只构建并联树，1：构建并联树及串联树
```

完成配置后，运行命令：run.py -c run.config开始运行。

## 4.结果文件说明

```
out/
├── 01_cds_format #格式化后的cds文件
├── 02_pep #cds对应pep文件
├── 03_hmm_out #hmmscan搜索结果
├── 04_OG #筛选的同源低/单拷贝基因OG
├── 05_seq #OG对应的序列
│   ├── 01_cds #cds序列
│   └── 02_pep #pep序列
├── 06_aln #序列比对结果
│   ├── 01_aln #序列比对结果
│   ├── 02_trim #序列修剪结果
│   └── 03_trim_rename #序列修剪后将基因id重命名为物种名
├── 07_tree #树
│   ├── 01_coatree #并联树
│   └── 02_contree #串联树
└── 08_result #结果，串联树，并联树
```
