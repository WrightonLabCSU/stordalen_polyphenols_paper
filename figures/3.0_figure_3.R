library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(readxl)
library(dplyr)
library(tidyr)
library(ggnewscale)
library(cowplot)

#####
# generating annotated bacterial and archaeal genome trees for figure 3
#####

# read in GTDB trees (from data_processing)
bac_tree = read.tree("5.4_gtdbtk.bac120.decorated.tree")
bac_tree
arc_tree = read.tree("5.3_gtdbtk.ar53.decorated.tree")
arc_tree

# read in taxonomy (from data_processing)
annotation = read.delim("5.2_97_classification_wOutgroup.txt",sep="\t",header = FALSE)
colnames(annotation)=c("MAG","tax")
annotation=annotation%>%separate(tax,into=c("d","p","c","o","f","g","s"),sep=";",remove=FALSE)

# read in phylum color codes
p_colors=read_excel("Supplementary_Table_3.xlsx",sheet = "color_codes")

# read in mag metaG data
mag_encoded_paths=read_excel("Supplementary_Table_2.xlsx",sheet="MAG_metaG_paths")
mag_encoded_paths= mag_encoded_paths %>%
  select(mag,path,presence)
mag_encoded_paths$pathway=mag_encoded_paths$path

# read in and create MAG metaT data
pp_metaT_mags = read_excel("Supplementary_Table_2.xlsx",sheet="MAG_metaT_polyphenol_genes")
info = read_excel("Supplementary_Table_2.xlsx",sheet="polyphenol_genes")
info= info %>%select(gene_id,gene_description,Transformation,Substrate,oxygen)
metaT_genes=inner_join(info,pp_metaT_mags,by=c("gene_id"="ID"))

mag_metaT_paths=inner_join(metaT_genes,mag_encoded_paths,by=c("MAG"="mag","Transformation"="path"))

# determine presence or absense of polyphenol expression in a site by each MAG
mag_activity=mag_metaT_paths%>%
  select(MAG,sample,geTMM)%>%
  mutate(site=substr(sample,1,1))%>%
  group_by(MAG,site)%>%
  summarise(sum=sum(geTMM))%>%
  mutate(pa=ifelse(sum>0,1,0))
# check how many MAGs in 1, 2, or 3 sites
mag_activity%>%ungroup()%>%group_by(MAG)%>%summarise(sum=sum(pa))%>%ungroup()%>%group_by(sum)%>%summarise(n=n())

# get number of pathways encoded per MAG
paths=read_excel("Supplementary_Table_2.xlsx", sheet="MAG_metaG_paths")
paths_summary = paths %>%
  select(mag,path,family,presence)%>%
  distinct()%>%
  group_by(mag,family)%>%
  summarise(count=sum(presence))

# make bacterial tree
bac_dat=as.data.frame(bac_tree$tip.label)
colnames(bac_dat)=c("MAG")
bac_dat=left_join(bac_dat,paths_summary,by=c("MAG"="mag"))
bac_dat=left_join(bac_dat,annotation,by="MAG")
bac_p_fill=left_join(bac_dat,p_colors,by=c("p"="phylum"))%>%select(p,code)%>%distinct()
bac_fill=as.data.frame(sort(bac_p_fill$p))
colnames(bac_fill)=c("p")
bac_fill=left_join(bac_fill,bac_p_fill,by="p")
mag_activity$site=factor(mag_activity$site,levels = c("P","S","E"))

#manually pull nodes for key genera
bac_parents=as_tibble(bac_tree)
nodeids=c(3111,2901,2862,2133,2319,3348,1534,3251,3239,3457,3484,3283,3276)
nodedf=data.frame(node=nodeids)
nodelab <- c("Terra","Palsa-89","Palsa-295","Palsa-504","Lustribacter","Palsa-883","UBA6159","Roseiarcus","Novosphing","Bog-1198","FEN-1197","Bog-931","Bradyrhiz")
poslist=c(1)
labdf <- data.frame(node=nodeids, label=nodelab, pos=poslist)

# undecorated tree
p=ggtree(bac_tree, layout="fan", branch.length = "branch.length", size=0.2,open.angle=50) 

# annotate tree
bac=p+
  geom_fruit(data=bac_dat,geom=geom_tile,
             mapping=aes(y=MAG,x=0,fill=p),width = 0.1,offset=0.02,show.legend=FALSE)+
  scale_fill_manual(values=bac_fill$code)+
  new_scale_fill() +
  geom_fruit(data=bac_dat,
             geom=geom_bar,
             mapping=aes(y=MAG,x=count),
             orientation="y",
             stat="identity", offset=0.02,grid.params = list(),axis.params = list(axis="x"))+
  new_scale_fill() + 
  geom_fruit(data=mag_activity,
             geom=geom_point,
             mapping=aes(y=MAG,x=site,size=pa/3,color=site),offset=0.01,grid.params = list(),axis.params = list(axis="x"))+
  geom_highlight(node=3111,fill="orange")+
  geom_highlight(node=2901,fill="red")+
  geom_highlight(node=2862,fill="orange")+
  geom_highlight(node=2133,fill="red")+
  geom_highlight(node=2319,fill="red")+
  geom_highlight(node=3348,fill="red")+
  geom_highlight(node=1534,fill="red")+
  geom_highlight(node=3251,fill="orange")+
  geom_highlight(node=3239,fill="red")+
  geom_highlight(node=3457,fill="red")+
  geom_highlight(node=3484,fill="red")+
  geom_highlight(node=3283,fill="red")+
  geom_highlight(node=3276,fill="red")+
  geom_highlight(node=2954,fill="yellow")+
  geom_highlight(node=2974,fill="yellow")+
  geom_highlight(node=2999,fill="yellow")+
  geom_highlight(node=3330,fill="yellow")+
  geom_highlight(node=3269,fill="yellow")+
  geom_highlight(node=2455,fill="yellow")+
  geom_highlight(node=3050,fill="yellow")+
  geom_highlight(node=3261,fill="yellow")+
  geom_highlight(node=3425,fill="yellow")+
  geom_highlight(node=2915,fill="yellow")+
  geom_highlight(node=2216,fill="yellow")+
  geom_highlight(node=2024,fill="yellow")+
  geom_highlight(node=3547,fill="yellow")+
  geom_highlight(node=3366,fill="yellow")

# make archaeal tree
arc_dat=as.data.frame(arc_tree$tip.label)
colnames(arc_dat)=c("MAG")
arc_dat=left_join(arc_dat,paths_summary,by=c("MAG"="mag"))
arc_dat=left_join(arc_dat,annotation,by="MAG")
arc_dat$novelty=ifelse(arc_dat$c=="c__","class",ifelse(arc_dat$o=="o__","order",ifelse(arc_dat$f=="f__","family",ifelse(arc_dat$g=="g__","genus",ifelse(arc_dat$s=="s__","species","mag")))))
arc_p_fill=left_join(arc_dat,p_colors,by=c("p"="phylum"))%>%select(p,code)%>%distinct()
arc_fill=as.data.frame(sort(arc_p_fill$p))
colnames(arc_fill)=c("p")
arc_fill=left_join(arc_fill,arc_p_fill,by="p")

# undecorated tree
a=ggtree(arc_tree, layout="fan", branch.length = "branch.length", size=0.2,open.angle=350) 
a

# annotate tree
arc_dat[nrow(arc_dat)+1,]<-c("20110800_E2D_12","lignans","36","d__Archaea;p__Thermoproteota;c__Bathyarchaeia;o__TCS64;f__TCS64;g__UBA8941;s_","d__Archaea","p__Thermoproteota","c__Bathyarchaeia","o__TCS64","f__TCS64","g__UBA8941","s_","species")
arc_dat$count=as.numeric(arc_dat$count)
arc=a+geom_fruit(data=arc_dat,geom=geom_tile,
                 mapping=aes(y=MAG,x=0,fill=p),width = 0.1,offset=0.02,show.legend=FALSE)+
  scale_fill_manual(values=arc_fill$code)+
  new_scale_fill() + 
  geom_fruit(data=arc_dat,
             geom=geom_bar,
             mapping=aes(y=MAG,x=count),
             orientation="y",
             stat="identity", offset=0.02,grid.params = list(),axis.params = list(axis="x"))+
  geom_fruit(data=mag_activity,
             geom=geom_point,
             mapping=aes(y=MAG,x=site,size=pa/3,color=site),offset=0.1,grid.params = list(),axis.params = list(axis="x"))+
  geom_highlight(node=90,fill="yellow")

plot_grid(bac,arc)
