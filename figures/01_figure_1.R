library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Hmisc)
library(corrplot)
library(cowplot)

set.seed(08211995)

# read in Data_S_1
data_s1_metadata=read_excel("Data_S1.xlsx",sheet="Samples_Metadata")
data_s1_enzyme_assays=read_excel("Data_S1.xlsx",sheet="Enzyme_Assays")
data_s1_fc=read_excel("Data_S1.xlsx",sheet="polyphenol_content")
data_s1_gas=read_excel("Data_S1.xlsx",sheet="porewater_gas")
data_s1_metaT=read_excel("Data_S1.xlsx",sheet="metaT_enzymes")

# merge sheets
data_s1=data_s1_metadata%>%
  left_join(.,data_s1_enzyme_assays,by="Sample")%>%
  left_join(.,data_s1_fc,by="Sample")%>%
  left_join(.,data_s1_gas,by="Sample")%>%
  left_join(.,data_s1_metaT,by="Sample")

palsa=data_s1%>%
  filter(Habitat=="Palsa")
bog=data_s1%>%
  filter(Habitat=="Bog")
fen=data_s1%>%
  filter(Habitat=="Fen")

# set colors
colors=c("#058000","#0001FF","#703C1B")

#####
# correlating data and generating plots for fig 1, fig SX
#####
# fe2 vs PPO assay
cor.test(data_s1$`fe2 (mM)`,data_s1$`PhenOx  (nmol activity /g/ hr)`,method="pearson") # p=1.231e-05
# just in palsa
cor.test(palsa$`fe2 (mM)`,palsa$`PhenOx  (nmol activity /g/ hr)`,method="pearson") # p=0.49
# just in bog
cor.test(bog$`fe2 (mM)`,bog$`PhenOx  (nmol activity /g/ hr)`,method="pearson") # p=0.1116
# just in fen
cor.test(fen$`fe2 (mM)`,fen$`PhenOx  (nmol activity /g/ hr)`,method="pearson") # p=0.03419
# plot
x1=data_s1%>%
  ggplot(aes(x=`fe2 (mM)`,y=`PhenOx  (nmol activity /g/ hr)`))+
  geom_smooth(method = "lm",color="#00aeef",fill="#00aeef")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
p1=palsa%>%
  ggplot(aes(x=`fe2 (mM)`,y=`PhenOx  (nmol activity /g/ hr)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#703C1B")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()
b1=bog%>%
  ggplot(aes(x=`fe2 (mM)`,y=`PhenOx  (nmol activity /g/ hr)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()
f1=fen%>%
  ggplot(aes(x=`fe2 (mM)`,y=`PhenOx  (nmol activity /g/ hr)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5, color="#0001FF")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()

# fe2 vs PPO metaT
cor.test(data_s1$`fe2 (mM)`,data_s1$PPO_metaT,method="pearson") # p=0.04976
fit=lm(fen$`fe2 (mM)`~fen$PPO_metaT)
summary(fit) 
# just in palsa
cor.test(palsa$`fe2 (mM)`,palsa$PPO_metaT,method="pearson") # p=0.04976
# just in bog
cor.test(bog$`fe2 (mM)`,bog$PPO_metaT,method="pearson") # p=0.01907
# just in fen
cor.test(fen$`fe2 (mM)`,fen$PPO_metaT,method="pearson") # p=0.2234

# plot
x2=data_s1%>%
  ggplot(aes(x=`fe2 (mM)`,y=PPO_metaT))+
  geom_smooth(method = "lm",color="#ed1c24",fill="#ed1c24")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
p2=palsa%>%
  ggplot(aes(x=`fe2 (mM)`,y=PPO_metaT))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#703C1B")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,4)
b2=bog%>%
  ggplot(aes(x=`fe2 (mM)`,y=PPO_metaT))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,4)
f2=fen%>%
  ggplot(aes(x=`fe2 (mM)`,y=PPO_metaT))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5, color="#0001FF")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,4)

## adjusting pvalues for correlations
# for overall sites by fe2
p=c(0.00001235,0.04976)
p.adjust(p,method = "BH") # 0.0000247 0.0497600
# within sites by fe2
p=c(0.49,0.9297,0.1116,0.01907,0.03419,0.2234)
p.adjust(p,method = "BH") #0.58800 0.92970 0.22320 0.10257 0.10257 0.33510

# PPO assay vs FC assay
cor.test(data_s1$`PhenOx  (nmol activity /g/ hr)`,data_s1$mgMeGal_dry,method="pearson") # p=0.2448
# just in palsa
cor.test(palsa$`PhenOx  (nmol activity /g/ hr)`,palsa$mgMeGal_dry,method="pearson") # p=0.7455
# just in bog
cor.test(bog$`PhenOx  (nmol activity /g/ hr)`,bog$mgMeGal_dry,method="pearson") # p=0.6868
# just in fen
cor.test(fen$`PhenOx  (nmol activity /g/ hr)`,fen$mgMeGal_dry,method="pearson") # p=0.2221
# plot
x3=data_s1%>%
  filter(Sample!="F_3_S",Sample!="F_3_D") %>%
  ggplot(aes(x=`PhenOx  (nmol activity /g/ hr)`,y=mgMeGal_dry))+
  geom_smooth(method = "lm",color="black")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
p3=palsa%>%
  ggplot(aes(x=`PhenOx  (nmol activity /g/ hr)`,y=mgMeGal_dry))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#703C1B")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.3)+
  ylim(0,450)
b3=bog%>%
  ggplot(aes(x=`PhenOx  (nmol activity /g/ hr)`,y=mgMeGal_dry))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.3)+
  ylim(0,450)
f3=fen%>%
  filter(Sample!="F_3_S",Sample!="F_3_D")%>%
  ggplot(aes(x=`PhenOx  (nmol activity /g/ hr)`,y=mgMeGal_dry))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5, color="#0001FF")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.3)+
  ylim(0,450)

# PPO metaT vs FTICRMS polyphenols
cor.test(data_s1$PPO_metaT,data_s1$percent_polyphenol_fticrms,method="pearson") # p=0.6973
# just in palsa
cor.test(palsa$PPO_metaT,palsa$percent_polyphenol_fticrms,method="pearson") # p=0.1949
# just in bog
cor.test(bog$PPO_metaT,bog$percent_polyphenol_fticrms,method="pearson") # p=0.2661
# just in fen
cor.test(fen$PPO_metaT,fen$percent_polyphenol_fticrms,method="pearson") # p=0.05331
# plot
x4=data_s1%>%
  ggplot(aes(x=PPO_metaT,y=percent_polyphenol_fticrms))+
  geom_smooth(method = "lm",color="black")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
p4=palsa%>%
  ggplot(aes(x=PPO_metaT,y=percent_polyphenol_fticrms))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#703C1B")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.8)+
  ylim(0,0.2)
b4=bog%>%
  ggplot(aes(x=PPO_metaT,y=percent_polyphenol_fticrms),)+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.8)+
  ylim(0,0.2)
f4=fen%>%
  filter(Sample!="F_3_S",Sample!="F_3_D")%>%
  ggplot(aes(x=PPO_metaT,y=percent_polyphenol_fticrms))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5, color="#0001FF")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.8)+
  ylim(0,0.2)

## adjusting pvalues for correlations
# for overall
p=c(0.2448,0.6973)
p.adjust(p,method = "BH") # 0.4896 0.6973
# within sites
p=c(0.7455,0.1949,0.6868,0.2661,0.2221,0.05331)
p.adjust(p,method = "BH") # 0.74550 0.39915 0.74550 0.39915 0.39915 0.31986


# FC assay vs RNA
cor.test(data_s1$mgMeGal_dry,data_s1$`RNA (ng/ul)`,method="pearson") # p=0.3797
# just in palsa
cor.test(palsa$mgMeGal_dry,palsa$`RNA (ng/ul)`,method="pearson") # p=0.4939
# just in bog
cor.test(bog$mgMeGal_dry,bog$`RNA (ng/ul)`,method="pearson") # p=0.561
# just in fen
cor.test(fen$mgMeGal_dry,fen$`RNA (ng/ul)`,method="pearson") # p=0.1038
# plot
x5=data_s1%>%
  ggplot(aes(x=mgMeGal_dry,y=`RNA (ng/ul)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
p5=palsa%>%
  ggplot(aes(x=mgMeGal_dry,y=`RNA (ng/ul)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#703C1B")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,450)+
  ylim(0,200)
b5=bog%>%
  ggplot(aes(x=mgMeGal_dry,y=`RNA (ng/ul)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,450)+
  ylim(0,200)
f5=fen%>%
  ggplot(aes(x=mgMeGal_dry,y=`RNA (ng/ul)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5, color="#0001FF")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,450)+
  ylim(0,200)

# FTICRMS polyphenol vs metaT richness
cor.test(data_s1$percent_polyphenol_fticrms,data_s1$metaT_richness,method="pearson") # p=0.03134
# just in palsa
cor.test(palsa$percent_polyphenol_fticrms,palsa$metaT_richness,method="pearson") # p=0.04525
# just in bog
cor.test(bog$percent_polyphenol_fticrms,bog$metaT_richness,method="pearson") # p=0.573
# just in fen
cor.test(fen$percent_polyphenol_fticrms,fen$metaT_richness,method="pearson") # p=6.161e-06
# plot
x6=data_s1%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=metaT_richness))+
  geom_smooth(method = "lm",color="black")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
p6=palsa%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=metaT_richness))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#703C1B")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,250)
b6=bog%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=metaT_richness))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,250)
f6=fen%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=metaT_richness))+
  geom_smooth(method = "lm",color="#00aeef",fill="#00aeef")+
  geom_point(size=2,alpha=0.5,color="#0001FF")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,250)

## adjusting pvalues for correlations
# for overall
p=c(0.3797,0.03134)
p.adjust(p,method = "BH") # 0.37970 0.06268
# within sites
p=c(0.4939,0.04525,0.561,0.573,0.1038,0.00000616)
p.adjust(p,method = "BH") # 0.57300000 0.13575000 0.57300000 0.57300000 0.20760000 0.00003696


# FC assay vs BG Assay
cor.test(data_s1$mgMeGal_dry,data_s1$`BG  (nmol activity /g/ hr)`,method="pearson") # p=0.7185
# just in palsa
cor.test(palsa$mgMeGal_dry,palsa$`BG  (nmol activity /g/ hr)`,method="pearson") # p=0.3029
# just in bog
cor.test(bog$mgMeGal_dry,bog$`BG  (nmol activity /g/ hr)`,method="pearson") # p=0.8735
# just in fen
cor.test(fen$mgMeGal_dry,fen$`BG  (nmol activity /g/ hr)`,method="pearson") # p=0.1328
# plot
x7=data_s1%>%
  ggplot(aes(x=mgMeGal_dry,y=`BG  (nmol activity /g/ hr)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
p7=palsa%>%
  ggplot(aes(x=mgMeGal_dry,y=`BG  (nmol activity /g/ hr)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#703C1B")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,450)+
  ylim(0,9000)
b7=bog%>%
  ggplot(aes(x=mgMeGal_dry,y=`BG  (nmol activity /g/ hr)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,450)+
  ylim(0,9000)
f7=fen%>%
  ggplot(aes(x=mgMeGal_dry,y=`BG  (nmol activity /g/ hr)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5, color="#0001FF")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,450)+
  ylim(0,9000)

# FTICRMS polyphenol vs metaT GH expression
cor.test(data_s1$percent_polyphenol_fticrms,data_s1$GH_metaT,method="pearson") # p=0.8633
# just in palsa
cor.test(palsa$percent_polyphenol_fticrms,palsa$GH_metaT,method="pearson") # p=0.1677
# just in bog
cor.test(bog$percent_polyphenol_fticrms,bog$GH_metaT,method="pearson") # p=0.2917
# just in fen
cor.test(fen$percent_polyphenol_fticrms,fen$GH_metaT,method="pearson") # p=0.0002735
# plot
x8=data_s1%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=GH_metaT))+
  geom_smooth(method = "lm",color="black")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
p8=palsa%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=GH_metaT))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#703C1B")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,200)
b8=bog%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=GH_metaT))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,200)
f8=fen%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=GH_metaT))+
  geom_smooth(method = "lm",color="#00aeef",fill="#00aeef")+
  geom_point(size=2,alpha=0.5, color="#0001FF")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,200)

## adjusting pvalues for correlations
# for overall
p=c(0.7185,0.8633)
p.adjust(p,method = "BH") # 0.8633 0.8633
# within sites
p=c(0.3029,0.1667,0.8735,0.2917,0.1328,0.0002735)
p.adjust(p,method = "BH") # 0.363480 0.333400 0.873500 0.363480 0.333400 0.001641


# FC assay vs % Carb
cor.test(data_s1$mgMeGal_dry,data_s1$percent_carbohydrate_fticrms,method="pearson") # p=0.03758
# just in palsa
cor.test(palsa$mgMeGal_dry,palsa$percent_carbohydrate_fticrms,method="pearson") # p=0.2482
# just in bog
cor.test(bog$mgMeGal_dry,bog$percent_carbohydrate_fticrms,method="pearson") # p=0.9102
# just in fen
cor.test(fen$mgMeGal_dry,fen$percent_carbohydrate_fticrms,method="pearson") # p=0.9384
# plot
x9=data_s1%>%
  ggplot(aes(x=mgMeGal_dry,y=percent_carbohydrate_fticrms))+
  geom_smooth(method = "lm",color="#ed1c24",fill="#ed1c24")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
p9=palsa%>%
  ggplot(aes(x=mgMeGal_dry,y=percent_carbohydrate_fticrms))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#703C1B")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,450)+
  ylim(0,0.09)
b9=bog%>%
  ggplot(aes(x=mgMeGal_dry,y=percent_carbohydrate_fticrms))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,450)+
  ylim(0,0.09)
f9=fen%>%
  ggplot(aes(x=mgMeGal_dry,y=percent_carbohydrate_fticrms))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5, color="#0001FF")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,450)+
  ylim(0,0.09)

# FTICRMS polyphenol vs FTICRMS percent carbohydrate
cor.test(data_s1$percent_polyphenol_fticrms,data_s1$percent_carbohydrate_fticrms,method="pearson") # p=0.001771
# just in palsa
cor.test(palsa$percent_polyphenol_fticrms,palsa$percent_carbohydrate_fticrms,method="pearson") # p=0.03138
# just in bog
cor.test(bog$percent_polyphenol_fticrms,bog$percent_carbohydrate_fticrms,method="pearson") # p=0.1107
# just in fen
cor.test(fen$percent_polyphenol_fticrms,fen$percent_carbohydrate_fticrms,method="pearson") # p=0.5141
# plot
x10=data_s1%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=percent_carbohydrate_fticrms))+
  geom_smooth(method = "lm",color="#ed1c24",fill="#ed1c24")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
p10=palsa%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=percent_carbohydrate_fticrms))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#703C1B")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,0.09)
b10=bog%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=percent_carbohydrate_fticrms))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,0.09)
f10=fen%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=percent_carbohydrate_fticrms))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5, color="#0001FF")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,0.09)

## adjusting pvalues for correlations
# for overall
p=c(0.001771,0.03758)
p.adjust(p,method = "BH") # 0.003542 0.037580
# within sites
p=c(0.03138,0.2482,0.1107,0.9102,0.5141,0.9384)
p.adjust(p,method = "BH") # 0.18828 0.49640 0.33210 0.93840 0.77115 0.93840

# FC assay vs porewater CO2
bog_fen=data_s1%>%
  filter(Habitat!="Palsa")%>%
  filter(`porewater CO2 (mM)`!="NA")
bog_fen$`porewater CO2 (mM)`=as.numeric(bog_fen$`porewater CO2 (mM)`)
bog_pw=bog_fen%>%filter(Habitat=="Bog")
fen_pw=bog_fen%>%filter(Habitat=="Fen")
cor.test(bog_fen$mgMeGal_dry,bog_fen$`porewater CO2 (mM)`,method="pearson") # p=0.02257
# just in bog
cor.test(bog_pw$mgMeGal_dry,bog_pw$`porewater CO2 (mM)`,method="pearson") # p=0.2529
# just in fen
cor.test(fen_pw$mgMeGal_dry,fen_pw$`porewater CO2 (mM)`,method="pearson") # p=0.6164
# plot
x11=bog_fen%>%
  ggplot(aes(x=mgMeGal_dry,y=`porewater CO2 (mM)`))+
  geom_smooth(method = "lm",color="#00aeef",fill="#00aeef")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
b11=bog_pw%>%
  ggplot(aes(x=mgMeGal_dry,y=`porewater CO2 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,310)+
  ylim(0,6)
f11=fen_pw%>%
  ggplot(aes(x=mgMeGal_dry,y=`porewater CO2 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5, color="#0001FF")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,310)+
  ylim(0,6)

# FTICRMS polyphenol vs porewater CO2
bog_fen=data_s1%>%
  filter(Habitat!="Palsa")%>%
  filter(`porewater CO2 (mM)`!="NA")
bog_fen$`porewater CO2 (mM)`=as.numeric(bog_fen$`porewater CO2 (mM)`)
bog_pw=bog_fen%>%filter(Habitat=="Bog")
fen_pw=bog_fen%>%filter(Habitat=="Fen")
# across sites
cor.test(bog_fen$percent_polyphenol_fticrms,bog_fen$`porewater CO2 (mM)`,method="pearson") # p=0.0312
# just in bog
cor.test(bog_pw$percent_polyphenol_fticrms,bog_pw$`porewater CO2 (mM)`,method="pearson") # p=0.6669
# just in fen
cor.test(fen_pw$percent_polyphenol_fticrms,fen_pw$`porewater CO2 (mM)`,method="pearson") # p=0.9865
# plot
x12=bog_fen%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=`porewater CO2 (mM)`))+
  geom_smooth(method = "lm",color="#00aeef",fill="#00aeef")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
b12=bog_pw%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=`porewater CO2 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,6)
f12=fen_pw%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=`porewater CO2 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5, color="#0001FF")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,6)

## adjusting pvalues for correlations
# for overall
p=c(0.02257,0.0312)
p.adjust(p,method = "BH") # 0.0312 0.0312
# within sites
p=c(0.2529,0.6669,0.6164,0.9865)
p.adjust(p,method = "BH") # 0.8892 0.8892 0.8892 0.9865

# FC assay vs porewater CH4
bog_fen=data_s1%>%
  filter(Habitat!="Palsa")%>%
  filter(`porewater CH4 (mM)`!="NA")
bog_fen$`porewater CH4 (mM)`=as.numeric(bog_fen$`porewater CH4 (mM)`)
bog_pw=bog_fen%>%filter(Habitat=="Bog")
fen_pw=bog_fen%>%filter(Habitat=="Fen")
cor.test(bog_fen$mgMeGal_dry,bog_fen$`porewater CH4 (mM)`,method="pearson") # p=0.8114
# just in bog
cor.test(bog_pw$mgMeGal_dry,bog_pw$`porewater CH4 (mM)`,method="pearson") # p=0.6065
# just in fen
cor.test(fen_pw$mgMeGal_dry,fen_pw$`porewater CH4 (mM)`,method="pearson") # p=0.2144
# plot
x13=bog_fen%>%
  ggplot(aes(x=mgMeGal_dry,y=`porewater CH4 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
b13=bog_pw%>%
  ggplot(aes(x=mgMeGal_dry,y=`porewater CH4 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,310)+
  ylim(0,0.3)
f13=fen_pw%>%
  ggplot(aes(x=mgMeGal_dry,y=`porewater CH4 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5, color="#0001FF")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,310)+
  ylim(0,0.3)

# FTICRMS polyphenol vs porewater CH4
bog_fen=data_s1%>%
  filter(Habitat!="Palsa")%>%
  filter(`porewater CH4 (mM)`!="NA")
bog_fen$`porewater CH4 (mM)`=as.numeric(bog_fen$`porewater CH4 (mM)`)
bog_pw=bog_fen%>%filter(Habitat=="Bog")
fen_pw=bog_fen%>%filter(Habitat=="Fen")
cor.test(bog_fen$percent_polyphenol_fticrms,bog_fen$`porewater CH4 (mM)`,method="pearson") # p=0.5091
# just in bog
cor.test(bog_pw$percent_polyphenol_fticrms,bog_pw$`porewater CH4 (mM)`,method="pearson") # p=0.7678
# just in fen
cor.test(fen_pw$percent_polyphenol_fticrms,fen_pw$`porewater CH4 (mM)`,method="pearson") # p=0.1927
# plot
x14=bog_fen%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=`porewater CH4 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
b14=bog_pw%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=`porewater CH4 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,0.3)
f14=fen_pw%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=`porewater CH4 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5, color="#0001FF")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,0.3)

## adjusting pvalues for correlations
# for overall
p=c(0.8114,0.5091)
p.adjust(p,method = "BH") # 0.8114 0.8114
# within sites
p=c(0.6065,0.7678,0.2144,0.1927)
p.adjust(p,method = "BH") # 0.7678 0.7678 0.4288 0.4288


# main text plots
plot_grid(x1,x2,
          x3,x4,
          x5,x6,
          x7,x8,
          x9,x10,
          x11,x12,
          x13,x14,ncol=2)
# SOM plots
plot_grid(p1,b1,f1,
          p2,b2,f2,
          p3,b3,f3,
          p4,b4,f4,
          p5,b5,f5,
          p6,b6,f6,
          p7,b7,f7,
          p8,b8,f8,
          p9,b9,f9,
          p10,b10,f10,
          "",b11,f11,
          "",b12,f12,
          "",b13,f13,
          "",b14,f14,ncol=6)
