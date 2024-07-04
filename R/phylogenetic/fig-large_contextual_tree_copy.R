
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



## DATA#################
### trees
# time scaled
all_phl.tree.time=read.beast("output/phylogenetic/dated_trees/all_wgsrate_lsd_CI.date.nexus")
# substitution scaled
all_phl.tree.subs=read.newick("output/phylogenetic/trees/all_rooted.nwk")

## metadata
all_phl.annot=read.csv("output/phylogenetic/concatenated_alignment/ph_metadata_merged_by_id_581.csv")

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
  gplot+aes(color=outbreak) +  layout_rectangular()+
  labs(caption="Time (years)")+
  scale_color_manual(values=c("black","darkred","blue"))+
  guides(colour = "none")+
  theme_tree2()+
  scale_x_continuous(breaks=seq(1890,2020, 10), minor_breaks=seq(1890, 2030, 5)) +
  geom_fruit(geom=geom_tile, mapping=aes(fill=outbreak), width=6, offset=0.05)+
  scale_fill_manual(values=c("transparent","darkred","blue"),
                    labels = c('','Romblon: 2022-23', 'Romblon: 2011-12'))+ theme(legend.title = element_blank())+new_scale_fill()+
  theme(panel.grid.major   = element_line(color="darkgrey", linewidth=.2,linetype="dotted"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 0.5)) + theme(legend.position = "top");tree.time



###substitution tree
# convert to substitutions:
all_phl.tree.subs$edge.length=all_phl.tree.subs$edge.length*11925
gplot2 <- ggtree(all_phl.tree.subs,ladderize = TRUE, size=0.1) %<+% all_phl.annot
# highlight the romblon cases 
tree.snps=
  gplot2+theme_tree2()+
  #geom_treescale(x=0, y=700, width=20, color='black', label="SNPs",linesize=1, offset=2)+
    labs(caption="Number of substitutions")+
aes(color=outbreak) +  layout_rectangular()+
  scale_color_manual(values=c("black","darkred","blue"))+
  guides(colour = "none")+
  geom_fruit(geom=geom_tile, mapping=aes(fill=outbreak), width=21)+
  scale_fill_manual(values=c("transparent","darkred","blue"),
                    labels = c('','Romblon: 2022-23', 'Romblon: 2011-12'))+ theme(legend.title = element_blank())+
  theme(panel.grid.major   = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  guides(fill="none");tree.snps

trees=tree.time|tree.snps 

trees+ plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 12))

tiff("output/figures/phylogenetic/S1_contextTrees_300dpi.tif", units = "in", width=11.69, height=8.27, res= 300, compression = "lzw")     
trees+ plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 12))
dev.off()

pdf("output/figures/phylogenetic/S1_contextTrees.pdf", width=11.69, height=8.27)     
trees+ plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 12))
dev.off()
