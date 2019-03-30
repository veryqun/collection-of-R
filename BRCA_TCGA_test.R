#########################
## 2019/3/29   4:20 PM ##
## 三阴性乳腺癌表行分析##
#########################
#三阴性乳腺癌是指癌组织免疫组织化学检查结果为雌激素受体（ER）、孕激素受体（PR）和原癌基因Her-2均为阴性的乳腺癌。
#这类乳腺癌占所有乳腺癌病理类型的10.0%～20.8%，具有特殊的生物学行为和临床病理特征，预后较其他类型差。

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



