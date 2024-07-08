
library(treeio)
library(phytools)
library(ggtree)
library(wesanderson)
library(ggplot2)  
library(ggpubfigs)
library(ggpubr)
library(ggtreeExtra)
library(RColorBrewer)
library(ape)
library(castor)
library(phangorn)
library(viridis)
library(ggnewscale)
library(patchwork)
library(ggstar)



## DATA#################
### trees
# time scaled
all_phl.tree.time=read.beast("output/phylogenetic/dated_trees/all_wgsrate_lsd_CI.date.nexus")
# substitution scaled
all_phl.tree.subs=read.newick("output/phylogenetic//trees/all_rooted.nwk")

## metadata
all_phl.annot=read.csv("output/phylogenetic/concatenated_alignment/ph_metadata_merged_by_id_581.csv")
all_phl.annot$outbreak[all_phl.annot$outbreak=="Other"]=NA

# transmission chain info
tc=read.csv("output/phylogenetic/transmission_trees/constree_epigen_prunedDT99_simlocs_chainids.csv")
tc$lineage_chain=paste0("chain",tc$lineage_chain)

# add tc info to metadata
all_phl.annot<- merge(all_phl.annot, tc, by = 'Sample.ID', all.x = TRUE)

# tree plots:  aesthetics, ladderize
gplot <- ggtree(all_phl.tree.time, mrsd='2023-03-01',ladderize = TRUE, size=0.1) %<+% all_phl.annot

# highlight romblon seq- new and old
Romblon.new=which(gplot$data$outbreak=="Romblon_new")
Romblon.old=which(gplot$data$outbreak=="Romblon_old")

# highlight the romblon cases 
tree.time=
  gplot +  layout_rectangular()+
  labs(caption="Time (years)")+
  theme_tree2()+
  scale_x_continuous(breaks=seq(1890,2020, 10), minor_breaks=seq(1890, 2030, 5)) +
  geom_fruit(
    geom=geom_star,
    mapping=aes(subset=(!is.na(outbreak)),fill=outbreak,  starshape=outbreak),
    position="identity",colour="black",
    starstroke=0.1, size=2
  )+
  scale_fill_manual(name = "Outbreak period",
                    labels =  c('Romblon: 2022-23', 'Romblon: 2011-12'),
values=c("darkred","blue"), na.translate=FALSE) +   
  scale_starshape_manual(name = "Outbreak period",
                         labels =  c('Romblon: 2022-23', 'Romblon: 2011-12'),
values=c(1,13), na.translate=FALSE)+
  guides(fill = guide_legend(order=1,nrow=2,override.aes = list(size = 5),title.position = "top"), starshape=guide_legend(order=1,nrow=2,title.position = "top"))+
  new_scale_fill()+
  geom_fruit(geom=geom_tile, mapping=aes(fill=lineage_chain), width=6)+
  scale_fill_manual(name = "Lineages",values=c("#008CF9", "#B80058","#EBAC23","#006E00", "#00BBAD"),labels = c('1','2','3','4','5'),na.translate=FALSE )+
  theme(panel.grid.major   = element_line(color="darkgrey", linewidth=.2,linetype="dotted"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 0.5)) + theme(legend.position = "top")+  guides(fill=guide_legend(nrow=2, title="Lineages",title.position = "top"));tree.time

###substitution tree
# convert to substitutions:
all_phl.tree.subs$edge.length=all_phl.tree.subs$edge.length*11925
gplot2 <- ggtree(all_phl.tree.subs,ladderize = TRUE, size=0.1) %<+% all_phl.annot
# highlight the romblon cases 
tree.snps=
  gplot2+theme_tree2()+
    labs(caption="Number of substitutions")+
 layout_rectangular()+
  geom_treescale(col="grey")+
  geom_fruit(
    geom=geom_star,
    mapping=aes(subset=(!is.na(outbreak)),fill=outbreak,  starshape=outbreak),
    position="identity",colour="black",
    starstroke=0.1, size=2
  )+  
  scale_fill_manual(name = "Outbreak period",
                        labels =  c('Romblon: 2022-23', 'Romblon: 2011-12'),
                        values=c("darkred","blue"), na.translate=FALSE) +
 guides(fill ="none", starshape="none")+
  new_scale_fill()+
  geom_fruit(geom=geom_tile, mapping=aes(fill=lineage_chain), width=20)+
  scale_fill_manual(name = "Lineages",values=c("#008CF9", "#B80058","#EBAC23","#006E00", "#00BBAD"),na.translate=FALSE )+
  theme(panel.grid.major   = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  guides(fill="none");tree.snps

trees=tree.time|tree.snps 

trees+ plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 12))

tiff("output/figures/phylogenetic/S1_contextTrees_alt1.tif", units = "in", width=11.69, height=8.27, res= 300, compression = "lzw")     
trees+ plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 12))
dev.off()

pdf("output/figures/phylogenetic/S1_contextTrees_alt1.pdf", width=11.69, height=8.27)     
trees+ plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 12))
dev.off()


##### 
# collapse related nodes for subtree
library(castor)
test=collapse_tree_at_resolution(all_phl.tree.time@phylo, resolution=4, rename_collapsed_nodes = T)$tree

# tree plots:  aesthetics, ladderize
zoom <- ggtree(test, mrsd='2023-03-01',ladderize = TRUE, size=0.1) %<+% all_phl.annot

# highlight the romblon cases 
zoom.tree=
  zoom +  layout_rectangular()+
  labs(caption="Time (years)")+
  theme_tree2()+
  scale_x_continuous(breaks=seq(1890,2020, 10), minor_breaks=seq(1890, 2030, 5)) +
  geom_fruit(
    geom=geom_star,
    mapping=aes(subset=(!is.na(outbreak)),fill=lineage_chain,  starshape=lineage_chain),
    position="identity",colour="black",
    starstroke=0.1, size=2
  )+
  scale_fill_manual(name = "Lineage",
                    labels =  c('2', '3','1,4,5'),
                    values=c("#B80058","#EBAC23","#008CF9"), na.translate=FALSE) +   
  scale_starshape_manual(name = "Lineage",
                         labels =  c('2', '3','1,4,5'),
                         values=c(1,1,1), na.translate=FALSE)+
  guides(fill = guide_legend(order=1,nrow=2,override.aes = list(size = 5),title.position = "top"), starshape=guide_legend(order=1,nrow=2,title.position = "top"))+
  theme(panel.grid.major   = element_line(color="darkgrey", linewidth=.2,linetype="dotted"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 0.5)) + theme(legend.position = c(0.2, 0.5));zoom.tree
