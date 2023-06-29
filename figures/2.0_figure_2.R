library(pheatmap)
library(viridis)
library(tidyr)
library(dplyr)

paths=read_excel("Supplementary_Table_2.xlsx",sheet="metaT_pathway_relabun")

####
# plotting Fig. 2 and Extended Data Fig. 5
####

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



ox=paths%>%select(-fig_2_number)%>%gather(-Transformation,key="Sample",value="abund")%>%separate(Transformation,into=c("path","ox"),sep=";")%>%group_by(Sample,ox)%>%summarise(sum=sum(abund))%>%mutate(site=substr(Sample,1,1))
anox=ox%>%filter(ox=="anoxic")%>%mutate(compartment=paste(substr(Sample,1,1),substr(Sample,3,3),sep=""))%>%
  mutate(sat=ifelse(compartment=="ED","sat",
                    ifelse(compartment=="EM","sat",
                           ifelse(compartment=="ES","sat",
                                  ifelse(compartment=="SD","sat",
                                         ifelse(compartment=="SM","sat","unsat"))))))

anoxic_plot=anox%>%
  ggplot()+
  geom_boxplot(aes(x=ordered(sat,levels=c("unsat","sat")),y=sum))+
  geom_jitter(aes(x=ordered(sat,levels=c("unsat","sat")),y=sum,color=site))+
  theme_classic()

# run ANOVA and store results
res_aov <- aov(sum ~ sat,data = anox)
# check for normality using the residuals
library(car)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals, id = FALSE)
# need the qqplot to be a straight line, normal distribution in histogram
shapiro.test(res_aov$residuals)
# p-value >0.05 indicates normality

#assumption 2:  are the variance equal?
leveneTest(sum ~ as.factor(site), data = anox)
# p-value >0.05 indicates equal variances

res_aov <- aov(sum ~ sat,data = anox)
summary(res_aov) # significant, p=0.000294
