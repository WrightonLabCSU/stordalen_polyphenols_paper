library(readxl)
library(dplyr)
library(tidyr)
library(writexl)

set.seed(08211995)

#read in mag-level metaT abundance data
mag_geTMM = read_excel("Data_S3.xlsx",sheet="metaT_MAGs")
# read in gtdb v2
tax = read_excel("Data_S3.xlsx",sheet="MAGs")
tax=tax%>%select(MAG,GTDB_v2_r207)
# join metaT abundance with taxonomy
mag_geTMM=left_join(mag_geTMM,tax,by="MAG")
# read in phlyum color codes
phylum=read_excel("Data_S3.xlsx",sheet="color_codes")

# read in MAG metaG abundance
abund = read_excel("Data_S3.xlsx",sheet="metaG_MAGs")
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
mag_encoded_paths=read_excel("Data_S2.xlsx",sheet="MAG_metaG_paths")
mag_encoded_paths= mag_encoded_paths %>%
  select(mag,path,presence)
#mag_encoded_paths$MAG_ID=paste(mag_encoded_paths$taxonomy,mag_encoded_paths$mag,sep=";")
mag_encoded_paths$pathway=mag_encoded_paths$path
#mag_encoded_paths= mag_encoded_paths %>%
#select(-family.x,-taxonomy)
  
pp_metaT_mags = read_excel("Data_S2.xlsx",sheet="MAG_metaT_polyphenol_genes")
info = read_excel("Data_S2.xlsx",sheet="polyphenol_genes")
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
# using Data_S2 tab polyphenol_genes, figure_2_rxn_number column, remove duplicate entries
# this will now be Data_S2 tab "metaT_pathway_relabun"
