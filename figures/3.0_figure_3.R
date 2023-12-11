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
p_colors=read_excel("Supplementary_Data_3.xlsx",sheet = "color_codes")

# read in mag metaG data
mag_encoded_paths=read_excel("Supplementary_Data_2.xlsx",sheet="MAG_metaG_paths")
mag_encoded_paths= mag_encoded_paths %>%
  select(mag,path,presence)
mag_encoded_paths$pathway=mag_encoded_paths$path

# read in and create MAG metaT data
pp_metaT_mags = read_excel("Supplementary_Data_2.xlsx",sheet="MAG_metaT_polyphenol_genes")
info = read_excel("Supplementary_Data_2.xlsx",sheet="polyphenol_genes")
info= info %>%select(gene_id,gene_description,Transformation,Substrate,oxygen)
metaT_genes=inner_join(info,pp_metaT_mags,by=c("gene_id"="ID"),relationship = "many-to-many")

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
paths=read_excel("Supplementary_Data_2.xlsx", sheet="MAG_metaG_paths")
paths_summary = paths %>%
  select(mag,path,family,presence)%>%
  distinct()%>%
  group_by(mag,family)%>%
  summarise(count=sum(presence))

# read in talent/dominance (from data_processing)
tal_dom=read.delim("4.2_talent_dominant.txt",sep="\t")

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

#pull nodes for key genera
bac_parents=as_tibble(bac_tree)
bac_tal_dom_nodes=bac_parents%>%mutate(label=gsub(":g__","=g__",label))%>%mutate(label=gsub("; g__","=g__",label))%>%separate(label,into=c("one","g"),sep="=")%>%filter(is.na(g)==0)%>%mutate(g=gsub("'","",g))%>%right_join(.,annotation,by=c("g"))%>%select(-MAG,-tax,-s,-one,-branch.length)%>%distinct()%>%mutate(genus=paste(d,p,c,o,f,g,sep=";"))%>%right_join(.,tal_dom,by="genus")%>%select(node,genus)%>%mutate(pos=1)
labdf=bac_tal_dom_nodes
colnames(labdf)=c("node","label","pos")

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
  geom_hilight(node=2024,fill="yellow")+ #g__Fen-1064, dominant
  geom_hilight(node=2036,fill="red")+ #g__Palsa-465, talent
  geom_hilight(node=2133,fill="orange")+ #g__Palsa-504
  geom_hilight(node=2197,fill="yellow")+ #g__Bog-473
  geom_hilight(node=2216,fill="yellow")+ #g__RAAP-2
  geom_hilight(node=2266,color="yellow")+ #g__Fen-455
  geom_hilight(node=2319,fill="red")+ #g__Lustribacter
  geom_hilight(node=2336,fill="yellow")+ #g__Palsibacter
  geom_hilight(node=2424,fill="yellow")+ #g__Palsa-1400
  geom_hilight(node=2455,fill="yellow")+ #g__UBA11358
  geom_hilight(node=2615,fill="yellow")+ #g__ELB16-189
  geom_hilight(node=2656,fill="yellow")+ #g__LD21
  geom_hilight(node=2686,fill="yellow")+ #g__Paludibacter
  geom_hilight(node=2862,fill="orange")+ #g__Palsa-295
  geom_hilight(node=2893,fill="yellow")+ #g__Palsa-187
  geom_hilight(node=2901,fill="red")+ #g__Palsa-89
  geom_hilight(node=2907,fill="yellow")+ #g__Bog-105
  geom_hilight(node=2915,fill="yellow")+ #g__Bog-159
  geom_hilight(node=2954,fill="yellow")+ #g__BOG-234
  geom_hilight(node=2974,fill="yellow")+ #g__Bog-209
  geom_hilight(node=2999,fill="yellow")+ #g__Sulfotelmatobacter
  geom_hilight(node=3050,fill="yellow")+ #g__TOLSYN
  geom_hilight(node=3104,fill="yellow")+ #g__CAIVDV01
  geom_hilight(node=3111,fill="orange")+ #g__Terracidiphilus
  geom_hilight(node=3230,fill="yellow")+ #g__2-02-FULL-66-22
  geom_hilight(node=1513,fill="yellow")+#f__Magnetospirillaceae;g__
  geom_hilight(node=3239,fill="red")+ #g__Novosphingobium
  geom_hilight(node=1533,fill="red")+ #g__BOG-934
  geom_hilight(node=1532,fill="red")+ #g__Rhodopila
  geom_hilight(node=1534,fill="red")+ #g__UBA6159
  geom_hilight(node=3251,fill="red")+ #g__Roseiarcus
  geom_hilight(node=3261,fill="yellow")+ #g__Methylocella
  geom_hilight(node=3264,fill="red")+ #g__Rhodoblastus
  geom_hilight(node=3269,fill="yellow")+ #g__Methylocystis
  geom_hilight(node=3276,fill="red")+ #g__Bradyrhizobium
  geom_hilight(node=3283,fill="red")+ #g__BOG-931
  geom_hilight(node=3330,fill="orange")+ #g__Acidocella
  geom_hilight(node=3348,fill="red")+ #g__Palsa-883
  geom_hilight(node=3361,fill="red")+ #g__AP-15
  geom_hilight(node=3366,fill="yellow")+ #g__BOG-933
  geom_hilight(node=3389,fill="yellow")+ #g__LMDS01
  geom_hilight(node=3413,fill="yellow")+ #g__2-12-FULL-64-23
  geom_hilight(node=3425,fill="yellow")+ #g__PALSA-1006
  geom_hilight(node=3457,fill="orange")+ #g__Bog-1198
  geom_hilight(node=3484,fill="orange")+ #g__FEN-1191
  geom_hilight(node=3586,fill="yellow")+ #g__PMG-095
  geom_hilight(node=3591,fill="yellow")+ #g__Fen-1137
  geom_hilight(node=468,color="yellow")+ #g__PALSA-693 actino
  geom_hilight(node=194,fill="yellow")+ #f__Ktedonobacteraceae;g__
  geom_hilight(node=1735,fill="yellow")+ #g__Fen-1087
  geom_hilight(node=3544,fill="yellow")+ #g__Fen-1135
  geom_hilight(node=3550,fill="yellow")+ #g__PNOF01
  geom_hilight(node=1049,color="yellow") #f__UBA7541;g__

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

#pull nodes for key genera
arc_parents=as_tibble(arc_tree)
arc_tal_dom_nodes=arc_parents%>%mutate(label=gsub(":g__","=g__",label))%>%mutate(label=gsub("; g__","=g__",label))%>%separate(label,into=c("one","g"),sep="=")%>%filter(is.na(g)==0)%>%mutate(g=gsub("'","",g))%>%right_join(.,annotation,by=c("g"))%>%select(-MAG,-tax,-s,-one,-branch.length)%>%distinct()%>%mutate(genus=paste(d,p,c,o,f,g,sep=";"))%>%right_join(.,tal_dom,by="genus")%>%select(node,genus)%>%mutate(pos=1)
labdf=arc_tal_dom_nodes
colnames(labdf)=c("node","label","pos")

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
  geom_highlight(node=90,fill="yellow")+ #g__Bog-38
  geom_highlight(node=40,fill="yellow")+ #g__Fen-7
  geom_highlight(node=95,fill="yellow")+ #g__Methanosarcina
  geom_highlight(node=98,fill="yellow") #g__Methanothrix

plot_grid(bac,arc)
