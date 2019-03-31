
# 实验设计
实验目的决定试验方法和途径。
**试验目的** :获取三阴性乳腺癌的正常组织和癌症组织的基因表达差异情况，比较三阴性乳腺癌中的基因表达变化情况。
**试验设计** :通过TCGA获取乳腺癌的RNA-seq表达数据，筛选出三阴性乳腺癌的样本，通过比较癌症和正常组织的表达差异。

# TCGA数据库简介
**一句话介绍**:TCGA数据库是一个由国家癌症研究所(National Cancer Institute)和美国人类基因组研究所(National Human Genome Research Institute)共同监督的一个项目。使用对患者样本的**高通量基因组测序**和分析技术来试图提供括基因表达谱，拷贝数变异分析，SNP基因分型，全基因组DNA甲基化分析，微RNA分析等信息。收录了33种癌症基因组测序数据。TCGA数据处理和整理比Oncoman和GEO困难一些。但是针对**肿瘤和癌症**所能提供的信息是很**完善和可靠**的。
TCGA和GEO存在的区别是，GEO存在各种研究领域和研究方向的NGC数据和分析。TCGA是专门针对肿瘤和癌症设立的。TCGA优势是丰富且规范的临床数据，以及针对每种癌型的大样本量。

# TCGA数据的获取
**背景介绍：**
三阴性乳腺癌是指癌组织免疫组织化学检查结果为雌激素受体（ER）、孕激素受体（PR）和原癌基因Her-2均为阴性的乳腺癌。这类乳腺癌占所有乳腺癌病理类型的10.0%～20.8%，具有特殊的生物学行为和临床病理特征，预后较其他类型差。--from:[百度百科:三阴性乳腺癌](https://baike.baidu.com/item/%E4%B8%89%E9%98%B4%E6%80%A7%E4%B9%B3%E8%85%BA%E7%99%8C/5716245)

**实验设计：**
三阴性乳腺癌的筛选标准是根据pheotype来确定的，表达量是通过RNAseq结果确定的，正常组织和癌症组织是通过病例号确定的。
1. TCGA项目的数据可以通过Genomic Data Commons Data Portal获取，即通过GDC来访问，访问地址：https://portal.gdc.cancer.gov/。
2. TCGA数据库公开免费，所以有许多针对TCGA数据进行整合的网站。还可以通过UCSC Xena进行下载：https://xena.ucsc.edu/public。试验数据的选择和目的息息相关，试验设计：
![GDC](https://img-blog.csdnimg.cn/20190331142234177.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80NDQ1MjE4Nw==,size_16,color_FFFFFF,t_70)
点击BRCA（乳腺癌）进入数据选择界面。
3. 选择`gene expression RNAseq`>`HTSeq - Counts`。注意不要用RPKM等经过了normlization的表达矩阵来分析。要使用Counts来进行差异分析，因为在差异分析时候会自动进行标准化。如果数据经过处理例如log2+1，则可以下载后逆运算转变回来。
4. 选择`phenotype`>`Phenotype`这里面有病人的病例信息等，可以通过统计筛选出三阴性乳腺癌的患者。
![BRCA](https://img-blog.csdnimg.cn/20190331142419617.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80NDQ1MjE4Nw==,size_16,color_FFFFFF,t_70)
5. 点击连接进去会看到详细的信息如下载地址、样品数、数据处理方法等。
![RNAseq界面详细信息](https://img-blog.csdnimg.cn/20190331151314133.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80NDQ1MjE4Nw==,size_16,color_FFFFFF,t_70)
# 数据预处理
**实验设计**：筛选出三阴性乳腺癌的患者ID，再筛选出同时有癌症组织样本和癌旁组织样本，计算初始Count值。
## 三阴性乳腺癌患者（TNBC）筛选
**实验设计**：在R语言中处理数据，选择breast_carcinoma_estrogen_receptor_status（ER）、PR、HER2受体一栏全为隐形（Negative）的患者。
原始的phenotype文件如下图，信息量巨大
![phenotype文件](https://img-blog.csdnimg.cn/20190331153244425.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80NDQ1MjE4Nw==,size_16,color_FFFFFF,t_70)
**R语言实现**：
### 读取文件
```shell  @R
phenotype_file <- read.table("TCGA-BRCA.GDC_phenotype.tsv",header = T, sep = '\t', quote = "")
phenotype_colnames <- as.data.frame(colnames(phenotype_file))
#读取文件并去读列名（penotype的类型）
table(phenotype_file$breast_carcinoma_estrogen_receptor_status)
#查看雌激素受体(ER)的表达状态。
table(phenotype_file$breast_carcinoma_progesterone_receptor_status)
#查看孕激素受体（PR）状态
table(phenotype_file$lab_proc_her2_neu_immunohistochemistry_receptor_status)
#查看原癌基因Her-2状态
```
![talbe](https://img-blog.csdnimg.cn/20190331154407709.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80NDQ1MjE4Nw==,size_16,color_FFFFFF,t_70)
### 列出三项指标的列表，方便筛选

```shell  @R
colnames_num <- grep("receptor_status",colnames(phenotype_file))
#在phenotype_file列中检索“receptor_status”,grep返回值是列名的列表中包含有“receptor_status”的列序号，因为三阴性乳腺癌的三项指标在这里都有phenotype_file字段

#>colnames_num
#[1] 20 25 67 77 87 93
phenotype_colnames <- colnames(phenotype_file)[colname_num]
#利用匹配的返回值取列名到phenotype_colnames中,phenotype_colnames中存储的是含有“receptor_status”字段的列名
eph <- phenotype_file[,colnames_num[1:3]]
#把三项指标筛选出来，建立一个新的表格。方便筛选全阴的行（ID号）
```
![eph查看](https://img-blog.csdnimg.cn/20190331155307921.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80NDQ1MjE4Nw==,size_16,color_FFFFFF,t_70)
函数学习：
```shell  @R
?apply
# b为：
# first second
# one       1      4
# two       2      5
# three     3      6
# apply(b,1,sum)
# 上面的指令代表对矩阵b进行行计算，分别对每一行进行求和。函数涉及了三个参数：
# 第一个参数是指要参与计算的矩阵；
# 第二个参数是指按行计算还是按列计算，1——表示按行计算，2——按列计算；
# 第三个参数是指具体的运算参数。
# 上述指令的返回结果为：
# one   two three 
# 5     7     9 
```
### TNBC的筛选
```shell  @R
tnbc_rownum <- apply(eph,1,function(x)sum(x=="Negative"))
#eph：记录三项指标类型的矩阵
#将eph按照列执行sum(x=="Negative"):统计各行中，eph列中有Negative的数量
#question
tnbc_rownum
tnbc_sample <- phenotype_file[tnbc_rownum == 3, 1]
#统计每一行中数值是3的（即全为“Negative”）并取第一列
tnbc_sample
save(tnbc_sample,file = 'tnbc_sample.Rdata')
#保存Rdata
```


