#########################
## 2019/3/29   4:20 PM ##
## 三阴性乳腺癌表行分析##
#########################
#三阴性乳腺癌是指癌组织免疫组织化学检查结果为雌激素受体（ER）、孕激素受体（PR）和原癌基因Her-2均为阴性的乳腺癌。
#这类乳腺癌占所有乳腺癌病理类型的10.0%～20.8%，具有特殊的生物学行为和临床病理特征，预后较其他类型差。
#请不要用RPKM等经过了normlization的表达矩阵来分析。


getwd()
setwd('/Users/chenchaoqun/SUPER_QUN/R_practice/TCGA_BRCA_GDC/')
getwd()
phenotype_file <- read.table("TCGA-BRCA.GDC_phenotype.tsv",header = T, sep = '\t', quote = "")
phenotype_colnames <- as.data.frame(colnames(phenotype_file))
#查看
table(phenotype_file$breast_carcinoma_estrogen_receptor_status)
#查看雌激素受体(ER)的表达状态。
table(phenotype_file$breast_carcinoma_progesterone_receptor_status)
#查看孕激素受体（PR）状态
table(phenotype_file$lab_proc_her2_neu_immunohistochemistry_receptor_status)
#查看原癌基因Her-2状态
colnames_num <- grep("receptor_status",colnames(phenotype_file))
#在phenotype_file列中检索“receptor_status”,grep返回值是列名的列表中包含有“receptor_status”的列序号
phenotype_colnames <- colnames(phenotype_file)[colname_num]
#利用匹配的返回值取列名到phenotype_colnames中
eph <- phenotype_file[,colnames_num[1:3]]
#把三项指标筛选出来，建立一个新的表格。方便筛选全阴的行（ID号）

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

####################
#筛选出三阴性乳腺癌#
####################
tnbc_rownum <- apply(eph,1,function(x)sum(x=="Negative"))
#eph：类型名矩阵
#将eph按照列执行sum(x=="Negative"):统计各行中，eph列中有Negative的数量
#question
tnbc_rownum
tnbc_sample <- phenotype_file[tnbc_rownum == 3, 1]
#统计每一行中数值是3的（即全为“Negative”）并取第一列
tnbc_sample
save(tnbc_sample,file = 'tnbc_sample.Rdata')
#保存Rdata


##表达矩阵得到以后不能直接使用

library(data.table)
#?data.table
#R中的data.table包提供了一个data.frame的高级版本，让你的程序做数据整型的运算速度大大的增加。
a <- fread("TCGA-BRCA.htseq_counts.tsv.gz",sep = '\t',header = T)
#读取表达矩阵文件
a <- as.data.frame(a)
#转换为数据框
a[1:4,1:4]
rownames(a) <- a[,1]
a <- a[,-1]
#去除第一列
genes <- row.names(a)
genes[1:10]

##接下来还原counts数，网站使用了log2(count+1)进行counts数转换，接下来进行还原

a <- a^2-1

##接下来是取样本，需要取118个三阴性乳腺癌的样本，并且是成对的样本，即既有癌组织又有癌旁组织。

tnbc_p <- substring(tnbc_sample,1,12)
#取tnbc的每个元素的第一个到第12个字符 
#question

all_p <- substring(colnames(a),1,12)
head(all_p)
#取表达矩阵的病理号前12位
table(all_p)
#这里如果病人号前12位相同，说明是即既有癌组织又有癌旁组织。
table(all_p) == 2
paired_p <- names(table(all_p)[table(all_p) == 2])
need_p <- intersect(tnbc_p,paired_p)
#need_p是118个三阴性乳腺癌的样本，并且是成对的样本。即既有癌组织又有癌旁组织
#intersect

exprSet <- a[,all_p %in% need_p]
#构建表达矩阵
table(apply(exprSet, 1, function(x){sum(x==0)<10}))
tmp <- apply(exprSet, 1, function(x){sum(x==0)<10})
#对exprSet进行按行执行函数function：
#如果表达矩阵是0则取1，在18个样本中，如果这个基因的表达矩阵超过10个样本都是0，则舍弃这个基因
#question
exprSet <- exprSet[tmp,]
save(exprSet,file = 'tnbc_paired_exprSet.Rdata')
#保存三阴性乳腺癌的双样本表达矩阵

##########################
## 2019/3/30   3:33 PM  ##
##  差异基因分析        ##
##########################

dim(exprSet)
class(exprSet)
str(exprSet)
#查看信息experSet

?substr
#提取或替换字符向量

#TCGA可以根据第14和第15位判断是癌组织还是癌旁组织。01表示癌症组织，11表示正常组织（?question癌旁组织）
colnames(exprSet)
#提取列名
substr(colnames(exprSet),14,15)
#提取列名字符串的第14和15位置字符
as.numeric(substr(colnames(exprSet),14,15))
#现在提取的字符变成数字 '11'和 11不等同
ifelse(as.numeric(substr(colnames(exprSet),14,15)) < 10,'tumor','normal')
#判断第14和15号位的数值<10 即等于01时候，返回tumor 否则返回normal
factor(ifelse(as.numeric(substr(colnames(exprSet),14,15)) < 10,'tumor','normal'))
#question 把它变成因子型
group_list=factor(ifelse(as.numeric(substr(colnames(exprSet),14,15)) < 10,'tumor','normal'))
table(group_list)
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

#差异基因分析

library(DESeq2)

#构建一个病例号和肿瘤分类的对应关系
colData <- data.frame(row.names = colnames(exprSet),group_list= group_list)
#question如果这里设置colData，则环境里面没有colData数据

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
#basemean是指该基因在所有样本中的均值
DESeq2_DEG <-  resOrdered
nrDEG=DESeq2_DEG[,c(2,6)]
#只提取logFC和padj
colnames(nrDEG)=c('log2FoldChange','pvalue')

##########################
## 2019/3/30   6:54 PM ##
##      注释基因       ##
##########################

library("clusterProfiler")
#加载注释包
library("org.Hs.eg.db")

colnames(DESeq2_DEG)[2] <- 'logFC'
colnames(DESeq2_DEG)[5] <- 'P.value'
nrDEG = DESeq2_DEG
#更改nrDEG的列名

logFC_cutoff<- with(nrDEG,mean(abs(logFC))+2*sd(abs(logFC)))
logFC_cutoff
mean(nrDEG$logFC)
#question 为什么返回结果是 NA
logFC_cutoff = 1.2
ifelse(nrDEG$P.value<0.05 & abs(nrDEG$logFC)>logFC_cutoff,ifelse(nrDEG$logFC>logFC_cutoff,"UP","DOWN"),'NOT')
#如果Padj<0.05并且logFC>1.2则是显著基因，否则是不显著基因，标记为NOT，如果显著上调标记为UP，下调标记为DOWN
nrDEG$change <- as.factor(ifelse(nrDEG$P.value<0.05 & abs(nrDEG$logFC)>logFC_cutoff,ifelse(nrDEG$logFC>logFC_cutoff,"UP","DOWN"),'NOT'))
#把标记变为factor存储为nrDEG的新列change中
table(nrDEG$change)
#查看显著上调下调，不显著的基因个数

keytypes(org.Hs.eg.db)
#查看org.Hs.eg.db中的ID分类和表示方式

library("stringr")
?str_sub

rownames(nrDEG)<- str_sub(rownames(nrDEG),1,15)
#使用把基因名取出来。去除小数点后面的数字
#question# substr(rownames(nrDEG),1,15)效果相同
#也可使用.分割列然后取前部分#question用哪个函数查一下#
nrDEG$ENSEMBL <- rownames(nrDEG)
#以行名新创建一列

?bitr
#生物id转换
dif <- bitr(rownames(nrDEG),fromType = 'ENSEMBL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
#question#  57.92% of input gene IDs are fail to map...
head(dif)
difs <- bitr(rownames(nrDEG),fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = org.Hs.eg.db)
#nrDEG$SYMBOL <- rownames(nrDEG)
nrDEG <- nrDEG[,-9]
?merge
a<- merge(nrDEG,df,by = 'ENSEMBL')
#question# Error in if (nx >= 2^31 || ny >= 2^31) stop("long vectors are not supported") : 
#missing value where TRUE/FALSE needed
#b<- nrDEG[nrDEG$change=='UP'|nrDEG$change=='DOWN',]





