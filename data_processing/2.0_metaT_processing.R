library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(edgeR)
library(stringr)

set.seed(08211995)

# count table output by htseq-count
data = read.delim("bbfiltered_95id_vEMERGE97_revstranded_gene_counts_rev.tsv",header = FALSE)
# file of gene lengths
lengths = read.delim("drep_genes_lengths.csv",header = FALSE,sep=",")

################
## step 1: filter out genes with less than 5 counts
################

# remove bottom rows that are summary stats
data =data[1:6604140,]
# add column names- you can check order by doing "ls *2016.cat.bbfiltered.mapped_vsEMERGEv1_95.FILTERED.NAME.SORTED.bam"
colnames(data)=c("gene","E1D","E1M","E1S","E2D","E2M","E2S","E3D","E3M","E3S","P1D","P1M","P1S","P2D","P2M","P2S","P3D","P3M","P3S","S1D","S1M","S1S","S2D","S2M","S2S","S3D","S3M","S3S")

# check sum of columns
colSums(data[,-1])
# add column of gene count sums
data$sum=rowSums(data[,2:28])
# check range of count sums
range(data$sum)

# removing genes with no counts
filtered_counts=data %>%
   filter(sum>0)
# now all should have at least 1 count
range(filtered_counts$sum)
row.names(filtered_counts)=filtered_counts[,1]
filtered_counts=filtered_counts[,-1]
# making counts <5 to 0
counts5 = as.data.frame(filtered_counts)
counts5[counts5==1|counts5==2|counts5==3|counts5==4] <-0
counts5$sum=rowSums(counts5[,1:27])
# check range of count sums
range(counts5$sum)
filtered_counts5=counts5 %>%
  filter(sum>0)
# check range of count sums - min should be 5
range(filtered_counts5$sum)

################
## step 2: calculate geTMM
################
filtered_counts5$gene=row.names(filtered_counts5)
# add column with gene lengths
x = inner_join(filtered_counts5,lengths,by=c("gene"="V1"))
rpk = ((x[,1:27]*10^3)/x$V2)
row.names(rpk)=x$gene
group <- c(rep("A",ncol(rpk)))
rpk.norm <- DGEList(counts=rpk,group=group)

# make table of library size
lib_size=c("E1D","E1M","E1S","E2D","E2M","E2S","E3D","E3M","E3S","P1D","P1M","P1S","P2D","P2M","P2S","P3D","P3M","P3S","S1D","S1M","S1S","S2D","S2M","S2S","S3D","S3M","S3S")
lib_size=as.data.frame(lib_size)
# manually adding
lib_size$size = c("2571306","2040287","1266088","4276426","4382832","1118279","2183480","2626471","195983","1481421","1462511","1395369","1970623","2085592","1671297"," 879084","1783621","1622658","11292225","12573321","2799265","11566839","5699147","931740","7995884","5918197","801370")
colnames(lib_size)=c("sample","lib.size")
lib_size$lib.size=as.numeric(lib_size$lib.size)
row.names(lib_size) = lib_size$sample

rpk.norm$samples$lib.size = lib_size[rownames(rpk.norm$samples),]$lib.size
rpk.norm$samples$lib.size
rpk.norm$samples

rpk.norm <- calcNormFactors(rpk.norm)
getmms <- cpm(rpk.norm)

getmms_df = as.data.frame(getmms)
range(getmms_df)
getmms_df$gene=row.names(getmms)

# add a column of MAG by splitting MAG name from gene names
getmms_df$MAG=getmms_df$gene
getmms_df$MAG=str_replace(getmms_df$MAG,"_c_","=")
getmms_df$MAG=str_replace(getmms_df$MAG,"_Ga_","=")
getmms_df$MAG=str_replace(getmms_df$MAG,"_ENA","=")
getmms_df$MAG=str_replace(getmms_df$MAG,"_k141","=")
getmms_df$MAG=str_replace(getmms_df$MAG,"_scaffold","=")
getmms_df$MAG=str_replace(getmms_df$MAG,"_NODE","=")
getmms_df=getmms_df%>%separate(MAG,into=c("MAG",NA),sep="=")

# write a table of per gene/sample geTMM values
write_xlsx(getmms_df, 'getmms_ge5_1X.xlsx')

# this data is stored in Zenodo repository: 10.5281/zenodo.7591900, tab "metaT_genes"
