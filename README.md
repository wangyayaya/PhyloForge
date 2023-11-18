# TreeTool使用说明

本程序可用于基于核基因组CDS序列、细胞器CDS序列、细胞器全基因组序列及SNP位点构建系统发育树

## 1.基于低拷贝核基因组CDS序列构建系统发育树

### 1.1单基因构树模式

#### 1.所需文件

该流程需要两个需要用户提供的文件：

1）需要做物种树的物种cds序列（建议每个基因只保留最长转录本），存放于一个文件夹下，物种cds文件建议以genus_species.cds (fa/fas/fasta)来命名，最终物种树中即以genus_species作为物种名。

2）与做物种树的物种对应的BUSCO核心基因集，该文件可以按如下步骤获取处理：

访问BUSCO数据库：[Index of /v5/data/lineages/ (ezlab.org)](https://busco-data.ezlab.org/v5/data/lineages/)或者https://busco-archive.ezlab.org/v3/ 下载与物种对应的数据集，例如被子植物一般下载有胚植物（Embryophyta odb10）核心基因集：

```shell
$ wget -c https://busco-archive.ezlab.org/v3/datasets/prerelease/embryophyta_odb10.tar.gz -O embryophyta_odb10.tar.gz
$ tar -zxvf embryophyta_odb10.tar.gz
$ cd embryophyta_odb10
$ cat hmms/*.hmm >embryophyta_odb10.hmm
$ hmmpress embryophyta_odb10.hmm
# 获取embryophyta_odb10.hmm路径，将其写入对应的配置文件中
$ realpath embryophyta_odb10.hmm
```

#### 2.依赖软件

1）hmmer软件，安装该软件后，主要使用hmmpress命令和hmmscan命令。目前大部分window版本的hmmer软件不能识别我们需要的hmm文件。**必需**

2）序列比对软件，可选择mafft或muscle，默认为mafft

3）序列修剪软件：trimal，**必需**

4）构树软件：可选择fasttree、iqtree或raxmal，默认为raxmal，最好是多线程版本，否则线程数设置不起作用

5）pal2nal，将蛋白比对序列转换为codon后做树，在配置文件中seq=codon时，需配置该软件

6）基因树合并软件：astral，**必需**

#### 3.设置软件配置文件

首先配置相关软件，可运行est -gsc >software.config生成软件配置文件

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
**#[lcn_opt]选项需要用到的软件必需要有对应软件，如[lcn_opt]中序列比对软件选择mafft，则必需配置mafft软件路径，muscle就无需配置**
**#如果选择使用codon序列构树，则软件需配置pal2nal软件**

运行treetool  -sc software.config命令配置软件

#### 4.运行参数配置说明

可通过treetool -grc >run.config获取配置文件

run.config参数说明:

```
[lcn_opt] #构建基于低拷贝/单拷贝核基因的树
in_path = D:\WANG\Species_Tree\cs\ #cds文件目录
out_path = D:\WANG\Species_Tree\Species_Tree\out\  #结果输出目录
thread = int(>=2) #线程数，默认10，最小2
orthodb = /path/to/buscoDB  #busco数据库路径，参考1.1
aln_software = muscle/mafft #序列比对软件，默认mafft
tree_software = iqtree/fasttree/raxmal(Multithreaded version) #构树软件，默认raxmal 
cover = 7 #每个OG中物种数量，当做树物种为100，该参数选择50时，即每个OG中物种覆盖率至少为50%。默认为物种数量，即物种覆盖率100%
copy_number = 1/2/3…… #基因拷贝数，1：单拷贝，2以上为低拷贝，默认1
mode = 0 #选择运行模式，默认0
#0：从头开始跑完全流程
#1：从头跑到筛选完OG（初步判断OG数目是否适当）
#2：单独运行筛选OG这一步骤（筛选出适当的OG数目）
#3：从筛选完OG后运行完后续全步骤
seq = cds/pep/codon/codon1/condon2 #选择用那种类型序列来构建物种树，默认cds.codon1:保留密码子第1，2位，codon2：保留密码子第三位
#注意，当选择使用codon做树时，并不是所有序列都能转换为codon，因此实际用来做树的OG数目可能会少于筛选到的数目，可使用ls out_path/06_aln/02_trim|wc -l查看实际做树的OG数目
coa_con = 0/1 #选择构建并联树或串联树，0：只构建并联树，1：构建并联树及串联树
```

完成配置后，运行命令：treetool -l run.config构建基于低拷贝核基因系统发育树；-s参数构建基于SNP位点的系统发育树；-o参数构建基于细胞器编码基因的系统发育树；-w。

Tips：

1. 基于低拷贝核基因组流程支持续跑，如需增删物种续跑，只需增删in_path中文件即可。如果需要修改某一物种的cds文件后续跑，且03_hmm_out下已生成对应的完整的文件，请删除03_hmm_out下对应的该物种结果文件再续跑，或将新的cds文件重命名。
2. 当seq选择condon时，有些不规范的序列无法转换为codon，会默认将该序列删除。实际建树的序列可能会少于筛选到的OG的数量。

#### 5.结果文件说明

```
out/
├── 01_cds_format #格式化后的cds文件
├── 02_pep #cds对应pep文件
├── 03_hmm_out #hmmscan搜索结果
├── 04_OG #筛选的同源低/单拷贝基因OG，当某一物种在该OG下有多个拷贝时，根据hmmscan结果，保留e值最小的一个
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
└── 08_result #结果，生成结果以串联/并联_seq_software.nwk格式命名
```

### **1.2多基因构树模式（这里还没完全整合进去，不确定有没有做的必要）**

在1.1单基因构树模式中，当copy_number参数设置大于1时，某一物种在某一OG下会有多个拷贝，单基因模式会根据hmmscan搜索结果，保留e值最小的一个基因。这里我们开发了多基因模式，即将有多个拷贝的基因全部保留，使用astral-pro(https://github.com/chaoszhang/A-pro)合并基因树，该模式只能构建并联树

所有设置同1.1，除了astral软件修改为astral-pro

```
[lcn_opt]
in_path = D:\WANG\Species_Tree\cs\ 
out_path = D:\WANG\Species_Tree\Species_Tree\out\
thread = int(>=2)
orthodb = /path/to/buscoDB
aln_software = muscle/mafft
tree_software = iqtree/fasttree/raxmal(Multithreaded version)
cover = 7
mode = 0
seq = cds/pep/codon
coa_con = 0
retain_multi_copy = T/F #当某一OG中某一物种基因有多个拷贝时，是否保留多基因建树，当选择T时，只能构建并联树，且astral软件需要astral-pro版本
```



## 2.基于细胞器（叶绿体等）基因组数据构建系统发育树

### 2.1 基于编码基因

#### 1.所需文件

该软件需要需要用户提供做物种树物种的细胞器cds序列，存放于一个目录下，物种cds文件建议以genus_species.cds (fa/fas/fasta)来命名，最终物种树中即以genus_species作为物种名。要求每个物种fasta文件中基因ID以基因名来命名（基因名不区分大小写），或者使用从ncbi下载而来的格式，如下列三种格式：

```
>psbA
ATGCCAACTATTAAACAACTTATTAGAAATACAAGACAGCCAATCAGAAATGTCACGAAATCCCCCGCTC
……
>psbA lcl|MH042531.1_cds_AXT17556.1_2 [gene=psbA] [protein=photosystem II protein D1] [protein_id=AXT17556.1] [location=complement(459..1520)] [gbkey=CDS]
ATGCCAACTATTAAACAACTTATTAGAAATACAAGACAGCCAATCAGAAATGTCACGAAATCCCCCGCTC
……
>lcl|MH042531.1_cds_AXT17556.1_2 [gene=psbA] [protein=photosystem II protein D1] [protein_id=AXT17556.1] [location=complement(459..1520)] [gbkey=CDS]
ATGACTGCAATTTTAGAGAGACGCGAAAGCGAAAGCCTATGGGGTCGCTTCTGTAACTGGATAACTAGCA
……
```

Tips：NCBI下载的数据应带有[gene=psbA]，以识别基因，否则以fasta文件中的基因ID作为基因名

#### 2.依赖软件

见1.2中2-6

#### 3.设置软件配置文件

同1.3

#### 4.运行参数配置说明

```
[organelle_opt]
in_path =  D:\WANG\Species_Tree\cs\  #cds文件目录
out_path =  D:\WANG\Species_Tree\out\  #输出结果目录
thread = 10
cover = 8
aln_software = mafft
tree_software = raxml
seq = codon
```



#### 5.结果文件说明

```
out/
├── 01_cds_format #格式化后的cds文件
├── 02_pep #cds对应的pep文件
├── 03_OG_all #根据基因名获得的所有OG
├── 04_OG #根据物种覆盖度筛选后的OG
#以下同1.5中的结果
├── 05_seq
├── 06_aln
│   ├── 01_aln
│   ├── 02_trim
│   └── 03_trim_rename
├── 07_tree
│   ├── 01_coatree
│   └── 02_contree
└── 08_result
```

### 2.2 基于细胞器全基因组序列

#### 1.所需文件

该软件需要需要用户提供做物种树物种的细胞器全基因组序列，物种全基因组序列文件建议以genus_species.fa(fas/fasta)来命名最终结果中即以genus_species作为物种名。

#### 2.依赖软件

见1.2中2-5及6

#### 3.设置软件配置文件

同1.3

#### 4.运行参数配置说明

```
[whole_genome_opt]
in_path = whole_genome_data
out_path = w_out
thread = 10
aln_software = mafft
tree_software = raal
```



#### 5.结果文件说明

```
out/
├── 01_sequence_format
│   └── seq.fa
├── 02_aln
│   ├── 01_aln
│   │   └── seq.aln
│   └── 02_trim
│       └── seq.trim
├── 03_tree
│   ├── RAxML_bestTree.WholeGenome
│   ├── RAxML_bipartitionsBranchLabels.WholeGenome
│   ├── RAxML_bipartitions.WholeGenome
│   ├── RAxML_bootstrap.WholeGenome
│   └── RAxML_info.WholeGenome
└── 04_result
    ├── WholeGenome_RAxML_bipartitionsBranchLabels.nwk
    └── WholeGenome_RAxML_bipartitions.nwk
```

## 3.基于SNP位点构建系统发育树

### 1.所需文件

过滤后的VCF格式文件

### 2.依赖软件

与配置文件对应的构树软件：treebest或者PhyML

### 3.设置软件配置文件

同1.3

### 4.运行参数配置说明

```
[snp_opt] #构建基于SNP位点的树
vcf_file = in.vcf
out_path = /path/to/output/
tree_software = treebst #构树软件，可选择treebest或者phyml
```

### 5.结果文件说明

```
out/
├── 转化格式后的fasta或phylip格式的文件
└── 构建的进化树
```

tips：所有测试数据均可在GitHub获取：[wangyayaya/EasySpeciesTree: CDS to species tree (github.com)](https://github.com/wangyayaya/EasySpeciesTree)
