library(readxl)
library(dplyr)
library(tidyr)
library(writexl)

set.seed(08211995)

#read in mag-level metaT abundance data
mag_geTMM = read_excel("Supplementary_Table_3.xlsx",sheet="metaT_MAGs")
# read in gtdb v2
tax = read_excel("Supplementary_Table_3.xlsx",sheet="MAGs")
tax=tax%>%select(MAG,GTDB_v2_r207)
# join metaT abundance with taxonomy
mag_geTMM=left_join(mag_geTMM,tax,by="MAG")
# read in phlyum color codes
phylum=read_excel("Supplementary_Table_3.xlsx",sheet="color_codes")

# read in MAG metaG abundance
abund = read_excel("Supplementary_Table_3.xlsx",sheet="metaG_MAGs")
# convert metaG abundance to tidy format
tidy_abund = abund %>%
  gather(-MAG,key="sample",value="tpm")
tidy_abund$site=substr(tidy_abund$sample,1,1)
tidy_abund$depth=substr(tidy_abund$sample,3,3)
# convert to avergae metaG abundance per MAG, per site/depth
ave_tidy_abund = tidy_abund %>%
  group_by(site,depth,MAG) %>%
  summarise(ave=mean(tpm))

# read in MAG encoded transformations
mag_encoded_paths=read_excel("Supplementary_Table_2.xlsx",sheet="MAG_metaG_paths")
mag_encoded_paths= mag_encoded_paths %>%
  select(mag,path,presence)
#mag_encoded_paths$MAG_ID=paste(mag_encoded_paths$taxonomy,mag_encoded_paths$mag,sep=";")
mag_encoded_paths$pathway=mag_encoded_paths$path
#mag_encoded_paths= mag_encoded_paths %>%
#select(-family.x,-taxonomy)

pp_metaT_mags = read_excel("Supplementary_Table_2.xlsx",sheet="MAG_metaT_polyphenol_genes")
info = read_excel("Supplementary_Table_2.xlsx",sheet="polyphenol_genes")
info= info %>%select(gene_id,gene_description,Transformation,Substrate,oxygen)
metaT_genes=inner_join(info,pp_metaT_mags,by=c("gene_id"="ID"))
tidy_metaT_genes = metaT_genes %>%
  select(Transformation,Substrate,oxygen,MAG,sample,geTMM,site,depth) %>%
  group_by(Transformation,Substrate,oxygen,MAG,geTMM,site,depth) %>%
  summarise(ave=mean(geTMM))
tidy_metaT_paths = tidy_metaT_genes %>%
  ungroup() %>%
  group_by(Transformation,Substrate,oxygen,MAG,site,depth) %>%
  summarise(sum=sum(ave))

#combine metaG potential with metaT expression
data1 = left_join(ave_tidy_abund,mag_encoded_paths,by=c("MAG"="mag"))
#add metaT 
data1$path=tolower(data1$path)
tidy_metaT_paths$Transformation=tolower(tidy_metaT_paths$Transformation)
data2 = left_join(data1,tidy_metaT_paths,by=c("MAG","site","depth","path"="Transformation"))
data3 = data2 %>% select(-path,-Substrate,-oxygen)
# add column about whether it's encoded in the MAG
data3$metaG=ifelse(data3$ave>0,ifelse(is.na(data3$presence)==0,"metaG","not"),"absent")
# if it's encoded in metaG and "on" in metaT, give geTMM otherwise say not on
data3$metaT=ifelse(is.na(data3$presence)==1,0,ifelse(is.na(data3$sum)==0,data3$sum,0))
# add column of "metaT", "metaG", or 0
data3$fill=ifelse(data3$metaT>0,1,ifelse(data3$metaG=="metaG",0.5,0))
data3$fill_code=ifelse(data3$metaT>0,ifelse(data3$metaG=="metaG","metaT","none"),ifelse(data3$metaG=="metaG","metaG","none"))
data3$compartment=paste(data3$site,data3$depth,sep="")
data3$sum_2= ifelse(data3$fill_code=="metaT",data3$sum,0)
data3$ave_2= ifelse(data3$fill_code=="metaG",data3$ave,0)
mag_data = data3 %>%
  ungroup()%>%
  select(compartment,MAG,ave_2,sum_2) %>%
  group_by(MAG,compartment) %>%
  summarise(metaG=mean(ave_2),metaT=sum(sum_2))
mag_data$fill=ifelse(mag_data$metaT>0,1,ifelse(mag_data$metaG>0,0.5,0))
mag_data = mag_data %>%
  ungroup() %>%
  select(-metaG,-metaT) %>%
  spread(key=compartment,value=fill)
# reorder columns
mag_data=mag_data[,c(1,7,6,5,10,9,8,4,3,2)]
mag_data=as.data.frame(mag_data)
row.names(mag_data)=mag_data[,1]
mag_data=mag_data[,-1]

# plotting sum pathway expression
mag_encoded_paths$path=tolower(mag_encoded_paths$path)
metaT_genes$Transformation=tolower(metaT_genes$Transformation)
mag_metaT_paths=inner_join(metaT_genes,mag_encoded_paths,by=c("MAG"="mag","Transformation"="path"))
ave_path_data = mag_metaT_paths %>%
  ungroup() %>%
  group_by(Transformation,Substrate,oxygen,sample) %>%
  summarise(sum=sum(geTMM)) %>%
  mutate(site=substr(sample,1,1),depth=substr(sample,3,3)) %>%
  ungroup() %>%
  group_by(site,depth,Substrate,Transformation) %>%
  summarise(ave=mean(sum),sd=sd(sum)) %>%
  mutate(compartment=paste(site,depth,sep=""))
ave_path_data$compartment=factor(ave_path_data$compartment,levels=c("PS","PM","PD","SS","SM","SD","ES","EM","ED"))

# getting pathway expression = making pathway expression on a scale of 0-1 per pathway, across samples
max=mag_metaT_paths %>%
  ungroup() %>%
  group_by(Transformation,Substrate,oxygen,sample) %>%
  summarise(sum=sum(geTMM)) %>%
  summarise(max=max(sum))
paths = mag_metaT_paths %>%
  ungroup() %>%
  group_by(Transformation,Substrate,oxygen,sample) %>%
  summarise(sum=sum(geTMM)) %>%
  left_join(.,max,by=c("Transformation","Substrate","oxygen")) %>%
  ungroup()%>%group_by(Transformation,Substrate,oxygen,sample) %>%
  summarise(rel=sum/max) %>%
  ungroup()%>%
  select(-Substrate,-oxygen) %>%
  spread(key=sample,value=rel)
paths[is.na(paths)] <- 0

# write table 
write_xlsx(paths,"pathway_table_relabun.xlsx")
# using Supplementary_Table_2 tab polyphenol_genes, figure_2_rxn_number column, remove duplicate entries
# this will now be Supplementary_Table_2 tab "metaT_pathway_relabun"


##########
# calculating talent
##########
g_pathways=mag_metaT_paths%>%left_join(.,tax,by="MAG")%>%separate(GTDB_v2_r207,into=c("d","p","c","o","f","g","s"),sep=";")%>%mutate(genus=paste(d,p,c,o,f,g,sep=";"))%>%filter(geTMM>0)%>%select(genus,site,pathway)%>%distinct()%>%group_by(genus,site)%>%summarise(n=n())

# figuring out the 95% threshold for pathways expressed (talent cut-off)
quantile(g_pathways$n, probs = 0.95)  
mean(g_pathways$n) #5.6
median(g_pathways$n) #4
sd(g_pathways$n)
density(g_pathways$n)

a=ggplot(g_pathways, aes(x = n, y = 1, fill = (0.5 - abs(0.5 - stat(ecdf)))<0.05)) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_discrete(name = "Tail probability <0.05")+
  geom_vline(xintercept = 15)+
  theme_classic()+xlab("Unique Pathways Expressed")+ylab("")

##########
# calculating dominance
##########
category=read_delim("4.1_pathway_categories.txt",sep="\t")

# removing false positive annotations
mag_metaT_paths2=mag_metaT_paths%>%filter(gene!="3300037104_19_Ga0395706_017518_1")%>%filter(gene!="3300037104_19_Ga0395706_013576_5")%>%filter(gene!="3300037104_19_Ga0395706_002833_4")%>%filter(gene!="3300037104_19_Ga0395706_002118_12")%>%filter(gene!="20110800_E1S_4_c_000000002448_2")%>%filter(gene!="20120700_S1D_59_c_000000099081_15")%>%filter(gene!="3300037104_19_Ga0395706_001630_7")%>%filter(gene!="PLGY01_ENA|PLGY01000051|PLGY01000051.1_69")%>%filter(gene!="PMEE01_ENA|PMEE01000053|PMEE01000053.1_14")%>%filter(gene!="PMNG01_ENA|PMNG01000155|PMNG01000155.1_37")%>%filter(gene!="PMNG01_ENA|PMNG01000199|PMNG01000199.1_11")%>%filter(gene!="3300037104_19_Ga0395706_002511_1")%>%filter(gene!="20120800_E2X_6_c_000000002298_6")%>%filter(gene!="20120800_E2X_6_c_000000002427_2")%>%filter(gene!="20150700_E25_9_c_000000001274_3")%>%filter(gene!="20150700_E25_9_c_000000001337_8")
total=mag_metaT_paths2%>%left_join(.,tax,by="MAG")%>%separate(GTDB_v2_r207,into=c("d","p","c","o","f","g","s"),sep=";")%>%mutate(genus=paste(d,p,c,o,f,g,sep=";"))%>%filter(geTMM>0)%>%select(genus,sample,site,depth,pathway,geTMM,Substrate)%>%distinct()%>%group_by(genus,sample,site,depth,pathway,Substrate)%>%summarise(sum=sum(geTMM))%>%left_join(.,category,by="pathway")%>%filter(is.na(group)==0)%>%ungroup()%>%distinct()%>%group_by(sample,site,depth,group)%>%summarise(total=sum(sum))

dom=mag_metaT_paths2%>%
  left_join(.,tax,by="MAG")%>%separate(GTDB_v2_r207,into=c("d","p","c","o","f","g","s"),sep=";")%>%mutate(genus=paste(d,p,c,o,f,g,sep=";"))%>%filter(geTMM>0)%>%select(genus,sample,site,depth,pathway,geTMM,Substrate)%>%distinct()%>%group_by(genus,sample,site,depth,pathway,Substrate)%>%summarise(sum=sum(geTMM))%>%left_join(.,category,by="pathway")%>%filter(is.na(group)==0)%>%ungroup()%>%group_by(genus,sample,site,depth,group)%>%summarise(group_sample_sum=sum(sum))%>%left_join(.,total,by=c("sample","site","depth","group"))%>%mutate(abun=group_sample_sum/total)%>%ungroup()%>%group_by(genus,site,depth,group)%>%summarise(ave_abun=mean(abun))
# figuring out the 95% threshold for expression (dominance cut-off)
quantile(dom$ave_abun, probs = 0.95)  
mean(dom$ave_abun) #5.6
median(dom$ave_abun) #4

b=ggplot(dom, aes(x = ave_abun*100, y = 1, fill = (0.5 - abs(0.5 - stat(ecdf)))<0.05)) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_discrete(name = "Tail probability <0.05")+
  geom_vline(xintercept = 10)+
  theme_classic()+xlab("Relative Contribution to Community Polymer, Monomer, or Phenolic/Benzoic Expression (%)")+ylab("")


# making som figure 
plot_grid(a,b,ncol = 1)


dom_pass=dom%>%filter(ave_abun>=0.10)%>%ungroup()%>%select(genus)%>%distinct()%>%mutate(dominant="y")
tal_pass=g_pathways%>%filter(n>=15)%>%ungroup()%>%select(genus)%>%distinct()%>%mutate(talent="y")
tal_dom=full_join(dom_pass,tal_pass,by="genus")
write.table(tal_dom,"4.2_talent_dominant.txt", sep="\t")






