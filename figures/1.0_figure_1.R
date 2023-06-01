library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
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

# add saturation code
data_s1=data_s1%>%
  mutate(sat=ifelse(Habitat=="Fen","sat",ifelse(Sample=="S1D"|Sample=="S2D"|Sample=="S3D"|Sample=="S1M"|Sample=="S2M"|Sample=="S3M","sat","unsat")))

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

# PO assay vs FC assay
cor.test(data_s1$`PhenOx  (nmol activity /g/ hr)`,data_s1$mgMeGal_dry,method="pearson") # p=0.2448
# just in palsa
cor.test(palsa$`PhenOx  (nmol activity /g/ hr)`,palsa$mgMeGal_dry,method="pearson") # p=0.7455
# just in bog
cor.test(bog$`PhenOx  (nmol activity /g/ hr)`,bog$mgMeGal_dry,method="pearson") # p=0.6868
# just in fen
cor.test(fen$`PhenOx  (nmol activity /g/ hr)`,fen$mgMeGal_dry,method="pearson") # p=0.2221
# plot
x1=data_s1%>%
  filter(Sample!="E3S",Sample!="E3D") %>%
  ggplot(aes(x=`PhenOx  (nmol activity /g/ hr)`,y=mgMeGal_dry))+
  geom_smooth(method = "lm",color="black")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
p1=palsa%>%
  ggplot(aes(x=`PhenOx  (nmol activity /g/ hr)`,y=mgMeGal_dry))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#703C1B")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.3)+
  ylim(0,450)
b1=bog%>%
  ggplot(aes(x=`PhenOx  (nmol activity /g/ hr)`,y=mgMeGal_dry))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.3)+
  ylim(0,450)
f1=fen%>%
  filter(Sample!="E3S",Sample!="E3D") %>%
  ggplot(aes(x=`PhenOx  (nmol activity /g/ hr)`,y=mgMeGal_dry))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5, color="#0001FF")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.3)+
  ylim(0,450)

# PO metaT vs FTICRMS polyphenols
cor.test(data_s1$PPO_metaT,data_s1$percent_polyphenol_fticrms,method="pearson") # p=0.6973
# just in palsa
cor.test(palsa$PPO_metaT,palsa$percent_polyphenol_fticrms,method="pearson") # p=0.1949
# just in bog
cor.test(bog$PPO_metaT,bog$percent_polyphenol_fticrms,method="pearson") # p=0.2661
# just in fen
cor.test(fen$PPO_metaT,fen$percent_polyphenol_fticrms,method="pearson") # p=0.05331
# plot
x2=data_s1%>%
  ggplot(aes(x=PPO_metaT,y=percent_polyphenol_fticrms))+
  geom_smooth(method = "lm",color="black")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
p2=palsa%>%
  ggplot(aes(x=PPO_metaT,y=percent_polyphenol_fticrms))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#703C1B")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.8)+
  ylim(0,0.2)
b2=bog%>%
  ggplot(aes(x=PPO_metaT,y=percent_polyphenol_fticrms),)+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.8)+
  ylim(0,0.2)
f2=fen%>%
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



# FC assay vs BG Assay
cor.test(data_s1$mgMeGal_dry,data_s1$`BG  (nmol activity /g/ hr)`,method="pearson") # p=0.7185
# just in palsa
cor.test(palsa$mgMeGal_dry,palsa$`BG  (nmol activity /g/ hr)`,method="pearson") # p=0.3029
# just in bog
cor.test(bog$mgMeGal_dry,bog$`BG  (nmol activity /g/ hr)`,method="pearson") # p=0.8735
# just in fen
cor.test(fen$mgMeGal_dry,fen$`BG  (nmol activity /g/ hr)`,method="pearson") # p=0.1328
# plot
x5=data_s1%>%
  ggplot(aes(x=mgMeGal_dry,y=`BG  (nmol activity /g/ hr)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
p5=palsa%>%
  ggplot(aes(x=mgMeGal_dry,y=`BG  (nmol activity /g/ hr)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#703C1B")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,450)+
  ylim(0,9000)
b5=bog%>%
  ggplot(aes(x=mgMeGal_dry,y=`BG  (nmol activity /g/ hr)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,450)+
  ylim(0,9000)
f5=fen%>%
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
x6=data_s1%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=GH_metaT))+
  geom_smooth(method = "lm",color="black")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
p6=palsa%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=GH_metaT))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#703C1B")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,200)
b6=bog%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=GH_metaT))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,200)
f6=fen%>%
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
x7=bog_fen%>%
  ggplot(aes(x=mgMeGal_dry,y=`porewater CO2 (mM)`))+
  geom_smooth(method = "lm",color="#00aeef",fill="#00aeef")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
b7=bog_pw%>%
  ggplot(aes(x=mgMeGal_dry,y=`porewater CO2 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,310)+
  ylim(0,6)
f7=fen_pw%>%
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
x8=bog_fen%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=`porewater CO2 (mM)`))+
  geom_smooth(method = "lm",color="#00aeef",fill="#00aeef")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
b8=bog_pw%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=`porewater CO2 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,6)
f8=fen_pw%>%
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
x9=bog_fen%>%
  ggplot(aes(x=mgMeGal_dry,y=`porewater CH4 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
b9=bog_pw%>%
  ggplot(aes(x=mgMeGal_dry,y=`porewater CH4 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,310)+
  ylim(0,0.3)
f9=fen_pw%>%
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
x10=bog_fen%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=`porewater CH4 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(aes(color=Habitat),size=2,alpha=0.5)+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()
b10=bog_pw%>%
  ggplot(aes(x=percent_polyphenol_fticrms,y=`porewater CH4 (mM)`))+
  geom_smooth(method = "lm",color="black")+
  geom_point(size=2,alpha=0.5,color="#058000")+
  guides(color="none",alpha="none",size="none")+
  theme_classic()+
  xlim(0,0.2)+
  ylim(0,0.3)
f10=fen_pw%>%
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


p=c(0.8114,0.5091,0.02257,0.0312,0.7185,0.8633)
p.adjust(p,method = "BH")


#####
## testing PO expression/activity between saturated samples
#####
x11=data_s1%>%
  ggplot(aes(x=sat,y=`PhenOx  (nmol activity /g/ hr)`))+
  geom_jitter(aes(color=Habitat),size=2,alpha=0.5)+
  geom_boxplot(fill="NA")+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()


#####
## Kruskal-Wallis ANOVA
#####
data_s1_PO=data_s1%>%
  filter(Sample!="E3S",Sample!="E3D")
# is there a significant difference in gas amount by treatment
kruskal.test(`PhenOx  (nmol activity /g/ hr)`~sat, data=data_s1)
# p-value = 0.06279
# if yes, run Dunn's test
dunnTest(value ~ treatment,
         data=data_subset,
         method="bh")

x12=data_s1%>%
  ggplot(aes(x=sat,y=PPO_metaT))+
  geom_jitter(aes(color=Habitat),size=2,alpha=0.5)+
  geom_boxplot(fill="NA")+
  guides(color="none",alpha="none",size="none")+
  scale_color_manual(values=colors)+
  theme_classic()


#####
## Kruskal-Wallis ANOVA
#####
# is there a significant difference in gas amount by treatment
kruskal.test(PPO_metaT~sat, data=data_s1)
# p-value = 0.005848
# if yes, run Dunn's test
dunnTest(PPO_metaT ~ sat,
         data=data_s1,
         method="bh")


# main text plots
plot_grid(x1,x2,
          x3,x4,
          x5,x6,
          x7,x8,
          x9,x10,
          x11,x12,ncol=2)
# SOM plots
plot_grid(p1,b1,f1,
          p2,b2,f2,
          p3,b3,f3,
          p4,b4,f4,
          p5,b5,f5,
          p6,b6,f6,
          "",b7,f7,
          "",b8,f8,
          "",b9,f9,
          "",b10,f10,ncol=6)
