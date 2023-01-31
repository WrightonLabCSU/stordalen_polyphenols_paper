library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)

set.seed(08211995)

# Read in Data_S2
data = read_excel("Data_S2.xlsx",sheet="polyphenol_genes")

# add column with total number of steps per pathway (sum of points per pathway)
data_steps= data %>% group_by(Transformation) %>% summarise(sum=sum(points))
#data=left_join(data,data_steps,by="header")
tidy_data = data %>%
  gather(-rxn,-gene_id,-gene_description,-Substrate,-Transformation,-subheader,-points,-oxygen,-figure_2_rxn_number,key="mag",value="count")
#make counts into present/absent
tidy_data$a = ifelse(tidy_data$count>0,1,0)

# determine pathway completion for each MAG
summarised_mag_paths = tidy_data %>%
  group_by(mag,Transformation,oxygen) %>%
  summarise(points_real=sum(points*a)) %>%
  left_join(.,data_steps,by="Transformation") %>%
  summarise(comp = points_real/sum)

#check range of completion- should be 0 to 1
range(summarised_mag_paths$comp)

#make back into a wide table
wide_mag_paths = summarised_mag_paths %>%
  spread(key=Transformation,value=as.numeric(comp))

# make column "count" that tracks the number of pathways >=50% complete per MAG
wide_mag_paths$count=rowSums(wide_mag_paths[,-1]>=0.5)
# filtering out MAGs with count = 0
filtered_drep97_mags = wide_mag_paths %>%
  filter(count>0)
pp_encoding = filtered_drep97_mags$mag

# get MAG taxonomy
tax=read_excel("Data_S3.xlsx",sheet="MAGs")
tax = tax %>%select(MAG, GTDB_v2_r207)
# make a presence absence data frame of pathways encoded by MAGs
pa = filtered_drep97_mags
pa2= as.data.frame(ifelse(pa[,2:101]>=0.5,1,0))
pa2$mag=filtered_drep97_mags$mag
pa3 = pa2 %>%
  gather(-mag,key="path",value="presence") %>%
  filter(presence>0)

info = data %>%select(Transformation, Substrate, oxygen)%>%distinct()

pa4 = left_join(pa3,info,by=c("path"="Transformation"))
pa4= pa4 %>%
  separate(Substrate,c("polyphenols","family","subfamily","compound"),sep=";")
pa4= left_join(pa4,tax,by=c("mag"="MAG"))

# write new XLSX table, this is Data_S2 tab "MAG_metaG_paths"
write_xlsx(pa4,"mag_pathwayge50per_wAA1.xlsx")
