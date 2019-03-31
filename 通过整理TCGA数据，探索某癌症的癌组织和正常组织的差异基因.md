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
**实验设计**：查找三项指标的列并列出
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
实验设计：查找指标的状态并统计；列出三阴性并筛选出来这些数据
```shell  @R
tnbc_rownum <- apply(eph,1,function(x)sum(x=="Negative"))
#eph：记录三项指标类型的矩阵
#将eph按照列执行sum(x=="Negative"):统计各行中，eph列中有Negative的数量
#通过查看阴性的数量，函数返回结果是一个数字向量。内容是0或1、2、3。注意数字的排列顺序是对应病人ID的行。
tnbc_rownum
tnbc_sample <- phenotype_file[tnbc_rownum == 3, 1]
#通过让 tnbc_rownum==3 判别，如果成立，会返回一个logical向量，logical向量是可以直接被矩阵引用的，只会读取TURE的行或着列
#统计每一行中数值是3的（即全为“Negative”）并取第一列即病人ID
#tnbc_sample：存储着三阴性乳腺癌的
tnbc_sample
save(tnbc_sample,file = 'tnbc_sample.Rdata')
#保存Rdata
```
![阴性指标统计](https://img-blog.csdnimg.cn/20190331161025151.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80NDQ1MjE4Nw==,size_16,color_FFFFFF,t_70)
## 基因表达矩阵的构建  
### 基因表达矩阵的读取和读取后格式修改
**实验设计**：转换data.frame ；行名修改；Count值还原
```shell  @R
library(data.table)
#?data.table
#R中的data.table包提供了一个data.frame的高级版本，让你的程序做数据整型的运算速度大大的增加。
a <- fread("TCGA-BRCA.htseq_counts.tsv.gz",sep = '\t',header = T)
#读取表达矩阵文件
a <- as.data.frame(a)
#转换为数据框，读取后，转换前的a的类型是"data.table" "data.frame"
a[1:4,1:4]
rownames(a) <- a[,1]
a <- a[,-1]
#去除第一列，行名修改。
genes <- row.names(a)
genes[1:10]
##接下来还原counts数，网站使用了log2(count+1)进行counts数转换，接下来进行还原
a <- a^2-1
#a存储的全是数值型向量，可以直接通过数学处理进行全表修改。
```
![RNAseq结果的读取与处理](https://img-blog.csdnimg.cn/20190331165126541.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80NDQ1MjE4Nw==,size_16,color_FFFFFF,t_70)
### 对表达矩阵进行筛选构建
**实验设计**：列出TNBC的ID（前半部分），所有病人ID，有双样品类型的seq的病人ID，取TNBC的ID和双样品类型的seq病人ID交集作为目的病人ID。最后用 all_p %in% need_p 构建logical向量，被原始表达矩阵引用，构建符合目的要求的表达矩阵
```shell  @R
##接下来是取样本，需要取118个三阴性乳腺癌的样本，并且是成对的样本，即既有癌组织又有癌旁组织。
tnbc_p <- substring(tnbc_sample,1,12)
#取tnbc的每个元素的第一个到第12个字符 
#tnbc是TNBC患者的ID前部分
all_p <- substring(colnames(a),1,12)
head(all_p)
#取表达矩阵的病理号前12位
table(all_p)
#这里如果病人号前12位相同，说明是即既有癌组织又有癌旁组织。
table(all_p) == 2
paired_p <- names(table(all_p)[table(all_p) == 2])
#去除所有有两个组织类型的病人ID前部分赋值给paired_id
need_p <- intersect(tnbc_p,paired_p)
#need_p是118个是成对的样本。即既有癌组织又有癌旁组织

exprSet <- a[,all_p %in% need_p]
#构建表达矩阵：三阴性乳腺癌的样本，并且是成对的样本，即既有癌组织又有癌旁组织
#这里通过 all_p %in% need_p 处理返回一个 logical向量，logical和第一个数据元素all_p一一对应，元素个数和顺序和all_P相同。向量可以在在data.frame中引用，只会读取 TRUE
```
![目的表达矩阵](https://img-blog.csdnimg.cn/20190331172815502.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80NDQ1MjE4Nw==,size_16,color_FFFFFF,t_70)
### 对构建的表达矩阵进行数据筛选
**实验设计**：利用apply函数统计条件符合情况，然后保存为logical向量，然后利用logical的被引用来筛选
```shell  @R
table(apply(exprSet, 1, function(x){sum(x==0)<10}))
tmp <- apply(exprSet, 1, function(x){sum(x==0)<10})
#对exprSet进行按行执行函数function：
#如果表达矩阵是0则取1，在18个样本中，如果这个基因的表达矩阵超过10个样本都是空值，则舍弃这个基因
#tmp是一个logical向量
exprSet <- exprSet[tmp,]
save(exprSet,file = 'tnbc_paired_exprSet.Rdata')
#保存三阴性乳腺癌的双样本表达矩阵
```

### 癌症组织和正常组织的区分和标记
实验设计：利用substr()函数截取患者ID的有效信息位置，用ifelse()函数进行判断，并用factor()函数生成因子型变量。最后赋值得到患者的样品类型
>Tips：TCGA可以根据第14和第15位判断是癌组织还是癌旁组织。01表示癌症组织，11表示正常组织
```shell  @R
#TCGA可以根据第14和第15位判断是癌组织还是癌旁组织。01表示癌症组织，11表示正常组织
colnames(exprSet)
#提取列名
substr(colnames(exprSet),14,15)
#提取列名字符串的第14和15位置字符
as.numeric(substr(colnames(exprSet),14,15))
#现在提取的字符变成数字 '11'和 11不等同
ifelse(as.numeric(substr(colnames(exprSet),14,15)) < 10,'tumor','normal')
#判断第14和15号位的数值<10 即等于01时候，返回tumor 否则返回normal
factor(ifelse(as.numeric(substr(colnames(exprSet),14,15)) < 10,'tumor','normal'))
#question 把它变成因子型，因子型是有顺序的。这里需要变换成因子型
group_list=factor(ifelse(as.numeric(substr(colnames(exprSet),14,15)) < 10,'tumor','normal'))
#group_list记录了表达矩阵中患者ID的对应的组织类型，和表达矩阵的顺序一样。
table(group_list)
```
### 表达矩阵的筛选（应该在上一步一起进行）
实验设计：修正小于0的（因为变换过程中可能会产生-1，会影响实验），检查并去除na值。
```shell  @R
exprSet <- ceiling(exprSet)
#?ceiling（）取整，不小于变量本身
#round()以四舍五入形式保存指定小数位数
#floor()不大于变量本身最大整数
#table(is.na(exprSet))
#pmax(),使得矩阵最小值是0？
save <- exprSet
exprSet[exprSet<0] <- 0
table(exprSet<0)
table(is.na(exprSet))
```
![group_list和colData](https://img-blog.csdnimg.cn/20190331195814303.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80NDQ1MjE4Nw==,size_16,color_FFFFFF,t_70)
# 患者的癌症组织和正常组织的基因差异表达分析
```shell  @R
library(DESeq2)
#构建一个病例号和肿瘤分类的对应关系
colData <- data.frame(row.names = colnames(exprSet),group_list= group_list)
#colData:是患者ID号和组织类型的对应关系

#构建DESeq()函数要求的表达式
dds <- DESeqDataSetFromMatrix(countData = exprSet,colData = colData,design = ~ group_list)
#countData = exprSet,指出DESeq的表达矩阵
#colData = colData,知名每个表达矩阵的分类，比如实验组&对照组，正常组织&癌症组织
#question
#design= ~group_list，因子型，指出不同组的区别，是有顺序的。
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~group_list)
dds <- DESeq(dds)
#进行差异基因分析
resultsNames(dds)
res <-  results(dds, contrast=c("group_list","tumor","normal"))
#用group_list来做引导文件，用tumor来比较normal组织

resOrdered <- res[order(res$padj),]
#把res差异分析文件通过padj来排序
head(resOrdered)
resOrdered=as.data.frame(resOrdered)
#把resOrdered变成数据框，
write.table(resOrdered, file="single cell DGE order.xls", sep="\t",quote=F)
#写出差异基因文件
```
![DEG文件](https://img-blog.csdnimg.cn/20190331200306640.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl80NDQ1MjE4Nw==,size_16,color_FFFFFF,t_70)

