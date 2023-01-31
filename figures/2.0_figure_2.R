library(pheatmap)
library(viridis)
library(tidyr)
library(dplyr)


paths=read_excel("Data_S2.xlsx",sheet="metaT_pathway_relabun")

# get the average relative expression per site/depth
ave_paths = paths%>%select(-fig_2_number)%>%
  gather(-Transformation,key="sample",value="sum")%>%
  mutate(compartment=paste(substr(sample,1,1),substr(sample,3,3),sep=""))%>%
  group_by(Transformation,compartment)%>%summarise(ave=mean(sum))%>%
  spread(key=compartment,value=ave)

ave_paths[is.na(ave_paths)] <- 0
ave_paths=as.data.frame(ave_paths)
row.names(ave_paths)=ave_paths[,1]
ave_paths=ave_paths[,-1]

#plotting z-score of average module relative expression
pheatmap(t(ave_paths),cluster_rows = FALSE,cluster_cols=TRUE,scale="column")
