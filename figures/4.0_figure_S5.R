library(readxl)
library(dplyr)
library(tidyr)
library(pheatmap)

set.seed(08211995)

# read in Data_S4
data = read_excel("Data_S4.xlsx",sheet="polyphenols_phenolics")

# find average value of each compound per habitat and depth compartment
ave_by_cmpart=data%>%select(-method)%>%gather(-name,key="sample",value="area")%>%
  mutate(compartment=paste(substr(sample,1,1),substr(sample,3,3),sep=""))%>%
  group_by(name,compartment)%>%summarise(ave=mean(area))%>%
  spread(key = name,value=ave)
ave_by_cmpart=as.data.frame(ave_by_cmpart)
row.names(ave_by_cmpart)=ave_by_cmpart[,1]
ave_by_cmpart=ave_by_cmpart[,-1]

# plot heatmap to be used in Fig. S5
pheatmap(ave_by_cmpart,scale="column",cluster_rows = FALSE)
