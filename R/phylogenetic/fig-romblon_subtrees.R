# subset contextual tree based on closest relatives of outbreak cases

rm(list=ls())  # This clears anything currently in the environment

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
all_phl.tree=read.beast("output/phylogenetic/pastml_analysis/all_mppa_province/named.tree_all_wgsrate_lsd_CI.date.nexus")

## metadata
all_phl.annot=read.csv("output/phylogenetic//concatenated_alignment/ph_metadata_merged_by_id_581.csv")

### alignments
all_phl.aln=read.fasta("output/phylogenetic/concatenated_alignment/ph_all_merged_by_id_581.fasta")

# transmission chain info
tc=read.csv("output/phylogenetic/transmission_trees/constree_epigen_prunedDT99_simlocs_chainids.csv")
tc$lineage_chain=paste0("chain",tc$lineage_chain)

# add tc info to metadata
all_phl.annot<- merge(all_phl.annot, tc, by = 'Sample.ID', all.x = TRUE)
write.csv(all_phl.annot,"output/phylogenetic/concatenated_alignment/ph_metadata_merged_by_id_581_withTC.csv", row.names = F)

# foundation tree plot:  aesthetics, ladderize
gplot <- ggtree(all_phl.tree, mrsd='2023-03-01',ladderize = TRUE, size=0.1) %<+% all_phl.annot

# highlight romblon seq- new and old
Romblon.new=which(gplot$data$outbreak=="Romblon_new")
Romblon.old=which(gplot$data$outbreak=="Romblon_old")
## note that 2 sequences are not present in the tree (were outliers in lsd so removed automatically):
all_phl.annot$Sample.ID[which(!all_phl.annot$Sample.ID %in% all_phl.tree@phylo$tip.label)]

###Get subtrees relevant to the Romblon outbreak

# slice large contextual tree by time
root_age = get_tree_span(all_phl.tree@phylo)$max_distance
slices=treeSlice(all_phl.tree@phylo, slice=root_age-3, trivial=T, prompt=FALSE)

# then identify the relevant subtrees
matching_subtrees <- list()

for (slice in 1:length(slices)){
 if(sum(gplot$data$label[Romblon.new] %in% slices[[slice]]$tip.label)>0){
   matching_subtrees[[slice]] <- slices[[slice]]
 }
}

matching_subtrees <- matching_subtrees[!sapply(matching_subtrees, is.null)]

# use these to find MRCAs in the main tree
mrca1=MRCA(all_phl.tree, matching_subtrees[[3]]$tip.label)
mrca2=MRCA(all_phl.tree, matching_subtrees[[1]]$tip.label)
mrca3=MRCA(all_phl.tree, matching_subtrees[[2]]$tip.label)

# and subset contextual tree to include Romblon cases plus at least 5 closest relatives
source("R/phylogenetic//tree_functions.R")
subset_tree1 <- subset_tree_with_min_additional_tips(all_phl.tree, mrca1, min_additional_tips = 5)
# tree2 enncompasses t1 tips:
subset_tree2 <- subset_tree_with_min_additional_tips(all_phl.tree, mrca2, min_additional_tips = 5)
# tree 3
subset_tree3<- subset_tree_with_min_additional_tips(all_phl.tree, mrca3, min_additional_tips = 5)

## Tree plots
# Cluster 1
cluster1=ggtree(subset_tree1, mrsd='2023-03-01',ladderize = TRUE, size=0.3) %<+% all_phl.annot
cluster1$data$ACRdates=format(as.Date(decimal2Date(as.numeric(cluster1$data$date))), "%d-%b-%y")

  cluster1.plot=cluster1+
  geom_tiplab(aes(label=id_case), size=2.5)+
 layout_rectangular()+
 # scale_color_manual(values=c("black",alpha("darkred",0.5),"blue"))+guides(fill = "none")+
  theme_tree2()+
    geom_fruit(geom=geom_tile, mapping=aes(fill=lineage_chain), width=1, offset=0.1)+
    scale_fill_manual(name = "Transmission chain",values=c("#008CF9", "#006E00", "#00BBAD","#EBAC23"),labels = c('1','4','5','3'),na.translate=FALSE )+
    guides(fill = FALSE)+
    new_scale_fill()+
    geom_fruit(geom=geom_tile, mapping=aes(fill=Province.y), width=1, offset=0.1)+
    scale_fill_manual(values=c("purple","pink",alpha("darkred",0.5)),
                      labels = c("Bulacan","Metropolitan Manila","Romblon"))+
    guides(fill = FALSE)+
     scale_x_continuous(breaks=seq(2000, 2023, 2), minor_breaks=seq(2000, 20))+
  theme(panel.grid.major   = element_line(color="darkgrey", size=.5,linetype="dotted"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())  +
    geom_nodelab(aes(subset = node == MRCA(subset_tree1, matching_subtrees[[3]]$tip.label), label=ACRdates) ,size=3,fontface=3, nudge_x = -1, nudge_y = 1)+
    geom_nodelab(aes(subset = node == MRCA(subset_tree1, matching_subtrees[[3]]$tip.label)+1, label=ACRdates), size=3,fontface=3, nudge_x = -1, nudge_y = 1)+
    geom_nodelab(aes(subset = node == MRCA(subset_tree1, matching_subtrees[[3]]$tip.label)+2, label=ACRdates), size=3,fontface=3,nudge_x = -1, nudge_y = 1)+
    geom_nodepoint(aes(subset = node == MRCA(subset_tree1, matching_subtrees[[3]]$tip.label)+1), col=alpha(alpha("darkred",0.5),0.5), size=3)+
    geom_nodepoint(aes(subset = node == MRCA(subset_tree1, matching_subtrees[[3]]$tip.label)+2), col=alpha("darkred",0.5), size=3)+
    geom_nodepoint(aes(subset = node == MRCA(subset_tree1, matching_subtrees[[3]]$tip.label)), col=alpha("darkred",0.5), size=3)+
    geom_nodepoint(aes(subset = node == MRCA(subset_tree1, MRCA(subset_tree1, subset_tree1@phylo$tip.label))), col="purple", size=3, pch=15 ); cluster1.plot
  

  
cluster2=ggtree(subset_tree2, mrsd='2023-03-01',ladderize = TRUE, size=0.3) %<+% all_phl.annot
cluster2$data$ACRdates=format(as.Date(decimal2Date(as.numeric(cluster2$data$date))), "%d-%b-%y")

cluster2.plot=
  #cluster2+
  scaleClade(cluster2, 83, .05)%>%
  ggtree::collapse(node=83,'max', fill=alpha("grey",0.4))+
  geom_tiplab(aes(label=id_case), size=2.5)+
  layout_rectangular()+
  theme_tree2()+
  geom_fruit(geom=geom_tile, mapping=aes(fill=lineage_chain), width=2.5, offset=0.1)+
  scale_fill_manual(values=c("#B80058","#EBAC23"),
                    labels = c('2','3'),na.translate=FALSE )+
  guides(fill = "none")+
  new_scale_fill()+
  geom_fruit(geom=geom_tile, mapping=aes(fill=Province.y), width=2.5, offset=0.1)+
  scale_fill_manual(values=c("deeppink","sienna3","purple","lightgreen","lavender","blue4","pink","tomato","steelblue1","darkseagreen",alpha("darkred",0.5),"cyan","coral","lightsalmon","gold1",alpha("darkred",0.5),"yellowgreen"),labels = sort(unique(cluster2$data$Province.y)))+
  guides(fill = "none")+
 scale_x_continuous(breaks=seq(1970, 2030, 5), minor_breaks=seq(2000, 2030, 1)) +
  theme(legend.title = element_blank())+
  theme(panel.grid.major   = element_line(color="darkgrey", size=.2,linetype="dotted"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 0.5))+
  geom_nodelab(aes(subset = node == 173, label=ACRdates), size=3,fontface=3, nudge_x = -6, nudge_y = 1)+
  geom_nodepoint(aes(subset = node == 173), col="black", size=3, pch=0); cluster2.plot

#+geom_strip("Z-22-121", "4B-23-04" ,"4B-23-03", barsize=2, color='red')

cluster2.plot=cluster2.plot+geom_text(x=2028, y=10,colour="#EBAC23",fontface="bold",label="| Transmission\n chain 3 cases", check_overlap = T, size=3);cluster2.plot
  


cluster3=ggtree(subset_tree3, mrsd='2023-03-01',ladderize = TRUE, size=0.3) %<+% all_phl.annot
cluster3$data$ACRdates=format(as.Date(decimal2Date(as.numeric(cluster3$data$date))), "%d-%b-%y")

cluster3.plot=cluster3+
  layout_rectangular()+
  geom_tiplab(aes(label=id_case), size=2.5)+
  theme_tree2()+
  geom_fruit(geom=geom_tile, mapping=aes(fill=lineage_chain), width=1, offset=0.1)+
  scale_fill_manual(name="Transmission chain",values="#EBAC23",
                    labels = '3',na.translate=FALSE )+
  guides(fill = FALSE)+
  new_scale_fill()+
  geom_fruit(geom=geom_tile, mapping=aes(fill=Province.y), width=1, offset=0.1)+
  scale_fill_manual(values=c("deeppink","blue4","lightsalmon",alpha("darkred",0.5)),
                    labels = sort(unique(cluster3$data$Province.y)) )+
  guides(fill = FALSE)+
  theme(panel.grid.major   = element_line(color="darkgrey", size=.5,linetype="dotted"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  scale_x_continuous(breaks=seq(2000, 2023, 2), minor_breaks=seq(2000, 20))+
  geom_nodelab(aes(subset = node == MRCA(subset_tree3, matching_subtrees[[2]]$tip.label), label=ACRdates), fontface=3,  nudge_x = -1.4, nudge_y = 0.6,size=3)+
  geom_nodepoint(aes(subset = node == MRCA(subset_tree3, matching_subtrees[[2]]$tip.label)), col=alpha("darkred",0.5), size=3)+
  geom_nodepoint(aes(subset = node == MRCA(subset_tree3, matching_subtrees[[2]]$tip.label)-1), col="deeppink", size=3, pch=15);cluster3.plot

# fake legend
# Create vectors for relevant provinces and transmission chains
relevant_provinces <- unique(c(cluster1$data$Province.y, cluster3$data$Province.y))
relevant_tc <- unique(c(cluster1$data$lineage_chain, cluster3$data$lineage_chain,cluster2$data$lineage_chain))

# Create expanded data frame with all combinations of province and lineage_chain
expanded_data <- expand.grid(province = relevant_provinces, lineage_chain = relevant_tc)


# Plot the data with ggplot
province.plot=ggplot(expanded_data, aes(x = province, fill =province)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(name="Province",values=c("deeppink","purple","blue4","pink","lightsalmon",alpha("darkred",0.5)),labels = c("Batangas","Bulacan","Laguna","Metropolitan Manila","Quezon","Romblon"), na.translate=T) +
  theme(legend.position = "top")+guides(fill=guide_legend(nrow=2,title.position="top"))

tc.plot=ggplot(expanded_data, aes(x = province, fill =lineage_chain)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(name = "Lineages",values=c("#008CF9", "#B80058","#EBAC23","#006E00", "#00BBAD"),labels = c('1','2','3','4','5'),na.translate=FALSE )+
  theme(legend.position = "top")+guides(fill=guide_legend(nrow=2,title.position="top"))

legend.province <- as_ggplot(get_legend(province.plot))
legend.tc <- as_ggplot(get_legend(tc.plot))
#legends=legend.province|legend.tc+  plot_layout(widths = c(3,0.5));legends
legends=legend.tc|legend.province+  plot_layout(widths = c(3,0.8));legends
source("R/phylogenetic/map_philippines_relevant_province.R")
source("R/phylogenetic/map_romblon_cluster1.R")
source("R/phylogenetic/map_romblon_cluster2.R")
source("R/phylogenetic/map_romblon_cluster3.R")
collect1=cluster1.plot+cluster1.map
collect3=cluster3.plot+cluster3.map
collect2=cluster2.plot+cluster2.map
plot.part1=((collect1/collect2/collect3)|philippines)& theme(panel.border = element_rect(colour = "black", fill=NA))
 

tiff("output/figures/phylogenetic/final_phylo_300dpi_panellab.tif", units = "in", width=11.69, height=8.27, res= 300, compression = "lzw")
legends/plot.part1+ plot_layout(heights = c(0.2,1)) +
  plot_annotation(tag_levels = list(c('', '', 'A1','A2','B1','B2','C1','C2','D')))&
  theme(plot.tag = element_text(size = 12))
dev.off()

pdf("output/figures/phylogenetic/final_phylo_panellab2.pdf", width=11.69, height=8.27)
legends/plot.part1+ plot_layout(heights = c(0.2,1)) +
  plot_annotation(tag_levels = list(c('', '', 'A1','A2','B1','B2','C1','C2','D')))&
  theme(plot.tag = element_text(size = 12))
dev.off()

jpeg("output/figures/phylogenetic/final_phylo_panellab2.jpg", quality=100)
legends/plot.part1+ plot_layout(heights = c(0.2,1)) +
  plot_annotation(tag_levels = list(c('', '', 'A1','A2','B1','B2','C1','C2','D')))&
  theme(plot.tag = element_text(size = 12))
dev.off()

