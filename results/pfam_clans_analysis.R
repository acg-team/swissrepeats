###############################
### NOTE: Code not working.
### In file src/disordered_repeats_annotate.py activate the idr_perc calculation. 
### therefore all files need to be recalculated...
### This could be done with the most recent version of swiss-prot/uniprot and the other DB and TRAL...
###
### idr_perc is calculated by the length of the disordered region divided by the length of the TR to which it overlaps.
### IDEA: - calculate this here in R from the data from the cluster (?). 
###       - get the pfam data from here: https://pfam.xfam.org/help#tabview=tab13
###       - connect the two and get some information out of it. (WHICH INFORMATION?????)



rm(list = ls(all = TRUE))
gc()
library("ggplot2")
library("plyr")
library("scales")
require(grid)
library(RColorBrewer)
library(stringr)
setwd("/home/matteo/polybox/MSc_ACLS/swissrepeat/results")
source("helpers.R")


sp_path = "data/swissprot_annotations.tsv"
sp_all = load_swissprot(sp_path)

# tr_consensus_disorder_path = "results/tr_annotations/tr_all_with_disorder.csv"
tr_consensus_disorder_path = "results/tr_idr_overlap/repeat_disorder_overlap_regions.csv"
tr_consensus_disorder = load_tr_annotations(tr_consensus_disorder_path)
tr_consensus_disorder_sp = merge(x = tr_consensus_disorder, y = sp_all, by = "ID", all.x = TRUE)

pfam_families=read.delim("/home/matteo/polybox/MSc_ACLS/swissrepeat/data/pfam_annotations.tab", header = TRUE, sep='\t')
pfam_clans=read.delim("/home/matteo/polybox/MSc_ACLS/swissrepeat/data/Pfam-A.clans.tsv", header = TRUE, sep='\t')

# Annotate protein with a clan
# no run, only to create file, takes long!
# add_clan <- function(proteinID){
#   pfam_families_str<-as.character(pfam_families[pfam_families$ID==proteinID, "Pfam_ID"])
#   pfam_families = unlist(strsplit(pfam_families_str, ";"))
#   clans<-c() 
#   for (family in pfam_families) {
#     if (family %in% pfam_clans$Pfam_ID) {
#       this_clan=as.character(pfam_clans[pfam_clans$Pfam_ID==family, "CL_ID"])
#       clans<-c(clans,this_clan)
#     }  
#   }
#   return(paste(clans, sep=";"))}
# 
# 
# if ("parallel" %in% rownames(installed.packages())){
#   # Multiprocessing by forking only if 'parallel' package is already installed
#   library(parallel)
#   ncores <- detectCores()-1
#   pfam_families$Pfam_clans=mclapply(pfam_families$ID, add_clan, mc.cores = ncores)
# } else{
#   pfam_families$Pfam_clans=lapply(pfam_families$ID, add_clan)
# }
# 
# pfam_families$Pfam_clans = vapply(pfam_families$Pfam_clans, paste, collapse = ", ", character(1L))
# write.table(pfam_families, file = "../data/pfam_annotations_with_clans.tab")

# load annotation file produced above:
pfam_families=read.delim("/home/matteo/polybox/MSc_ACLS/swissrepeat/data/pfam_annotations_with_clans.tab", header = TRUE, sep="")

############ Plotting ###############
tr_consensus_disorder_sp = merge(x = tr_consensus_disorder_sp, y = pfam_families, by = "ID", all.x = TRUE)
#tr_consensus_disorder_sp$Pfam_clans=lapply(as.character(tr_consensus_disorder_sp$ID), add_clan)
disorder_type_bin <- with(tr_consensus_disorder_sp, cut(idr_perc, breaks = c(-0.1,0, 0.1, 0.3, 0.5, 1.0)))
tr_consensus_disorder_sp$disorder_group = disorder_type_bin

sorted_clan_counts=sort(table(unlist(tr_consensus_disorder_sp[tr_consensus_disorder_sp$l_type=='micro'&tr_consensus_disorder_sp$idr_perc>=0.8,"Pfam_clans"])), decreasing=TRUE)
sorted_clan_counts_df = data.frame(x=names(sorted_clan_counts[1:30]), y=sorted_clan_counts[1:30], row.names=NULL) 
sorted_clan_counts_df = merge(sorted_clan_counts_df,pfam_clans[!duplicated(pfam_clans[,c('CL_Name','CL_ID')]),c('CL_Name','CL_ID')], by.x = "x", by.y = "CL_ID",all.x = TRUE)

sorted_clan_counts_df=na.omit(sorted_clan_counts_df[ order(-sorted_clan_counts_df[,"y"]), ])
ggplot(sorted_clan_counts_df[sorted_clan_counts_df$x!="",],aes(x=CL_Name, y=y.Freq)) +
  geom_bar(stat='identity') + scale_x_discrete(limits = sorted_clan_counts_df$CL_Name) +
  labs(title="Pfam clans, micro repeats", x="Pfam clan", y="Number of repeats") +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1))

