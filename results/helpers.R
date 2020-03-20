# install.packages("plyr")
# install.packages("pals") # nice color gradients
# install.packages("ggplot2")
# install.packages("scales")
# install.packages("grid")
# install.packages("RColorBrewer")
# install.packages("tidyverse")
# install.packages("ggrepel")
# install.packages("xtable")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")
# install.packages("doParallel")
# install.packages("foreach")
# install.packages("taxize")
# install.packages("UpSetR")
# install.packages("reshape2")

require("plyr")
require("pals") # nice color gradients
require("ggplot2")
require("scales")
require("grid")
require("RColorBrewer")
require("tidyverse")
require("ggrepel")
require("xtable")
require("Biostrings")
require("doParallel")
require("foreach")
require("taxize")
require("UpSetR")
require("reshape2")


# In the local config file, *local_base_path* is defined. *local_base_path* is the base path where 
# all results files that need to be read in are located.
# In the local config file, *local_path_separator* is defined. *local_path_separator* is most likely
# '/' on a Unix system and '\' on a Windows system.
source("local_config.R")

# Set colors
cols1.4 <- c("#2D882D", "#AA3939", "#AA7939", "#29506D", "#5C7881")
cols1.4.bright <- c("#87CC87", "#FFAAAA", "#FFDAAA", "#708DA4", "#B2BEC1")
cols2.4 <- c("#2D882D", "#AA3939", "#AA7939", "#29506D")

cols1 <- c("#AA3939", "#AA7939", "#29506D", "#2D882D")  #http://paletton.com/#uid=7000I0kllllaFw0g0qFqFg0w0aF
cols1.bright <- c("#FFAAAA", "#FFD0AA", "#669999", "#88CC88")
cols2 <- c("#FFAAAA", "#FFDBAA", "#718EA4", "#88CC88") 
cols3 <- c("#801515", "#805215", "#123652", "#116611")
cols4 <- c("#550000", "#553100", "#042037", "#004400")
morecolors <- c("#ffbfc8", "#ffaa00", "#9979f2", "#00f281", "#ced9a3", "#d96c6c", "#6dcc00", "#2d98b3", "#998273", "#008077", "#736d1d", "#4d5766", "#592816", "#45264d", "#40f2ff", "#ff0000", "#3d3df2", "#f20000", "#6cd9a6", "#d9368d", "#ccbe00", "#a60016", "#8c6969", "#7f4400", "#004d73", "#33663a", "#590024", "#00294d", "#ff0044", "#b6eef2", "#3d9df2", "#e55c00", "#d9a66c", "#d936ce", "#b386b0", "#2c00a6", "#46468c", "#731d6d", "#007300", "#665533", "#000c59", "#0d332b")
morecolors2 <- c("#ff0000", "#ffaa00", "#ffff00", "#00aaff", "#ff00aa", "#00f2f2", "#e50000", "#e5e600", "#d90000", "#d9d900", "#006cd9", "#00cc88", "#bfbf00", "#00bfbf", "#b2b300", "#009999", "#7f1500", "#7f4000", "#80006a", "#007339", "#006073", "#001373", "#4c0026", "#331100", "#332200", "#2b0033", "#ff8040", "#ff40ff", "#793df2", "#6cd936", "#b22d43", "#73641d", "#1a330d", "#0d2033", "#ff8095", "#73bfe6", "#bf8f60", "#465e8c", "#5e468c", "#ffffbf", "#ffbfca", "#ace6c9", "#e6ace6", "#b39586", "#8c8169", "#8c697b", "#394d43", "#403530")
# Order taken from Uversky's paper: http://www.tandfonline.com/doi/full/10.4161/idp.24684
# Not mentioned in Uversky's paper: "B", "O", "U", "Z", "X". These guys might need to fit in with the rest (if possible, as some of them represent multiple aa.)
aa_order_promoting_to_disorder_promoting = c("C", "W", "I", "Y", "F", "L", "H", "V", "N", "M", "R", "T", "D", "G", "A", "K", "Q", "S", "E", "P", "B", "O", "U", "Z", "X")

# Not in operator
'%!in%' <- function(x,y)!('%in%'(x,y))

# detect available cores:
numcores <- detectCores()-1


load_swissprot <- function(path, tr_all){
  path = paste(local_base_path, path, sep=local_path_separator)
  sp_all = read.csv(path, sep="\t", header=TRUE, quote="", stringsAsFactors = FALSE)
  sp_all = plyr::rename(sp_all, c("Entry"="ID", "Cross.reference..DisProt."="DisProt_ID", 
                            "Cross.reference..MobiDB."="MobiDB_ID", 
                            "Taxonomic.lineage..SUPERKINGDOM."="Superkingdom",
                            "Taxonomic.lineage..KINGDOM."="Kingdom",
                            "Taxonomic.lineage..CLASS."="Class",
                            "Taxonomic.lineage..ORDER."="Order",
                            "Taxonomic.lineage..FAMILY."="Family",
                            "Taxonomic.lineage..SPECIES."="Species",
                            "Taxonomic.lineage.IDs..SPECIES."="Species_ID",
                            "Cross.reference..OrthoDB."="OrthoDB_ID",
                            "Cross.reference..OMA."="OMA_ID",
                            "Protein.names"="protein_name",
                            "Virus.hosts"="virus_hosts"))
  
  # Text-based check for some plastids. Note that a few cases are labeled as both,
  # and this misses rare orgenelles (tonoplast, amyloplastic, anaplastic, apicoplastic, chromoplastic, plastidial)
  sp_all %<>% mutate(is_chloroplastic =  grepl("chloroplastic", protein_name, T) & Superkingdom != "Viruses",
                     is_mitochondrial =  grepl("mitochondrial", protein_name, T) & Superkingdom != "Viruses") %>%
    mutate(origin=if_else(is_mitochondrial & !is_chloroplastic, paste0("Mitochondrial (", if_else(Kingdom != "", Kingdom, "unknown"), ")"),
                   if_else(is_chloroplastic & !is_mitochondrial, paste0("Chloroplastic (", if_else(Kingdom != "", Kingdom, "unknown"), ")"),
                   # if_else(is_chloroplastic, if_else(is_mitochondrial, paste0("Chloroplastic/Mitochondrial (", if_else(Kingdom != "", Kingdom, "unknown"), ")"), ),
                   if_else(is_chloroplastic & is_mitochondrial, paste0("Chloroplastic/Mitochondiral (", if_else(Kingdom != "", Kingdom, "unknown"), ")"),
                   if_else(Superkingdom == "Eukaryota", paste0(if_else(Kingdom != "", Kingdom, "unknown"), " (", Superkingdom, ")"),
                   if_else(Superkingdom == "Viruses", "Virus",
                   Superkingdom))))))

  # Add boolean to sp_all: has_tr. Later, add: has TR of specific type. 
  if(! missing(tr_all)){
    sp_proteins_w_trs = unique(tr_all$ID)
    sp_proteins_w_homorep = unique(tr_all[tr_all$l_effective==1, 'ID'])
    sp_proteins_w_microsats = unique(tr_all[tr_all$l_effective<=3, 'ID'])
    sp_proteins_w_short_trs = unique(tr_all[tr_all$l_effective>3 & tr_all$l_effective<15, 'ID'])
    sp_proteins_w_domain_trs = unique(tr_all[tr_all$l_effective>=15, 'ID'])
    
    sp_all$has_tr = sp_all$ID %in% sp_proteins_w_trs
    sp_all$has_homo_tr = sp_all$ID %in% sp_proteins_w_homorep
    sp_all$has_micro_tr = sp_all$ID %in% sp_proteins_w_microsats
    sp_all$has_short_tr = sp_all$ID %in% sp_proteins_w_short_trs
    sp_all$has_domain_tr = sp_all$ID %in% sp_proteins_w_domain_trs  
  }
  
  return(sp_all)
}

load_tr_annotations <- function(path) {
  path = paste(local_base_path, path, sep=local_path_separator)
  tr_all = read.csv(path, header = TRUE, quote="")
  tr_all = subset(tr_all, pvalue < 0.01)
  tr_all$total_repeat_length = (tr_all$n_effective * tr_all$l_effective)
  tr_all$center = tr_all$begin + (tr_all$total_repeat_length - 1)/2
  tr_all$l_type = ifelse(tr_all$l_effective ==1, "homo", 
                         ifelse(tr_all$l_effective >1 & tr_all$l_effective <= 3, "micro", 
                                ifelse(tr_all$l_effective < 15, "small", 
                                       "domain")))
  tr_all$fraction_disordered_chars = tr_all$disordered_overlap / (tr_all$l_effective * tr_all$n_effective)
  return(tr_all)
}

load_overlap_annotations <- function(path) {
  path = paste(local_base_path, path, sep=local_path_separator)
  tr_all = read.csv(path, header = TRUE, quote="")
  tr_all = subset(tr_all, pvalue < 0.01)
  tr_all$total_repeat_length = (tr_all$n_effective * tr_all$l_effective)
  tr_all$center = tr_all$begin + (tr_all$total_repeat_length - 1)/2
  tr_all$l_type = ifelse(tr_all$l_effective <= 3, "micro", ifelse(tr_all$l_effective < 15, "small", "domain"))
  return(tr_all)
}
load_disorder_annotations <- function(path){
  path = paste(local_base_path, path, sep=local_path_separator)
  discoor_all = read.csv(path, header = TRUE, quote="")
  discoor_all = plyr::rename(discoor_all, c("uniprotID"="ID")) 
  return(discoor_all)
}

load_homorepeat_data <- function(path, aa_ignore, set){
  path = paste(local_base_path, path, sep=local_path_separator)
  data = read.csv(path, header = TRUE, quote="")
  data$set = as.factor(set)
  data$count_rounded = round(data$count)
  data$log10_count = log10(data$count)
  data$log10_count_rounded = log10(data$count_rounded)
  data$repeat_region_length = data$n * data$count
  data$aa = factor(data$aa, levels=aa_order_promoting_to_disorder_promoting)
  if(! missing(aa_ignore)){
    data = subset(data, !(aa %in% aa_ignore))
  }
  n_chars_total = sum(data[data$type=="empirical","repeat_region_length"])
  data$frequency = data$count/n_chars_total
  return(data)
}

load_expected_homorepeat_frequencies <- function(path, aa_ignore, set){
  path = paste(local_base_path, path, sep=local_path_separator)
  data = read.csv(path, header = TRUE, quote="")
  data$set = as.factor(set)
  data$aa = factor(data$aa, levels=aa_order_promoting_to_disorder_promoting)
  if(! missing(aa_ignore)){
    data = subset(data, !(aa %in% aa_ignore))
  }
  return(data)  
}


load_amino_acid_counts <- function(path, type, superkingdom, aa_ignore){
  path = paste(local_base_path, path, sep=local_path_separator)
  data = read.csv(path, header = TRUE, sep=",")
  # Transpose the data.frame
  aa = colnames(data)
  count = as.numeric(as.vector(data[1,]))
  data = data.frame(aa, count)
  
  # Ignore aas in aa_ignore
  data = subset(data, ! (aa %in% aa_ignore))
  
  # Add more info
  data$type = type
  data$Superkingdom = superkingdom
  data$frequency = data$count/sum(data$count)
  data$aa = factor(data$aa, levels=aa_order_promoting_to_disorder_promoting)
  return(data) 
}

compare_amino_acid_counts <- function(data_order, data_disorder){
  data_order$relative_frequency_ordered_divided_by_disordered = data_order$frequency / data_disorder$frequency
  data_order$log10relative_frequency_ordered_divided_by_disordered = log10(data_order$relative_frequency_ordered_divided_by_disordered)
  data_order$total_frequency = ((data_order$frequency) * sum(data_order$count) +  (data_disorder$frequency) * sum(data_disorder$count))/(sum(data_order$count) + sum(data_disorder$count))
  
  data_disorder$relative_frequency_ordered_divided_by_disordered = data_order$relative_frequency_ordered_divided_by_disordered
  data_disorder$log10relative_frequency_ordered_divided_by_disordered = data_order$log10relative_frequency_ordered_divided_by_disordered
  data_disorder$total_frequency = data_order$total_frequency
 
  return(rbind(data_order, data_disorder))
}

load_pfam_map <- function(path_to_pfam_map){
  # If pfam_mapping doesn't exist download it. Else load it.
  if(!file.exists(path_to_pfam_map)){
    res <- tryCatch(download.file("ftp://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/pdb_pfam_mapping.txt", 
                                  destfile = paste0(local_base_path, local_path_separator, "data/pdb_pfam_mapping.txt")),
                    error=function(e) 1)}
  if(file.exists(path_to_pfam_map)) {
    pfam_map <- read.csv(paste0(local_base_path, local_path_separator, "data/pdb_pfam_mapping.txt"), sep="\t", header = TRUE) 
  }
  # split the PFAM_ACC column by "." (has to be escaped with \\)
  split_cols <- colsplit(as.character(pfam_map$PFAM_ACC), "\\.",names = c("PFAM_ACC", "PFAM_ACC_suf"))
  # replace the original PFAM_ACC column with the new PFAM_ACC column and add the suffix as new column
  pfam_map$PFAM_ACC <- split_cols$PFAM_ACC
  pfam_map <- cbind(pfam_map, split_cols$PFAM_ACC_suf)
  names(pfam_map)[names(pfam_map) == "split_cols$PFAM_ACC_suf"] <- "PFAM_ACC_suf"
  return(pfam_map)
}

load_sp_protein_sequences <- function(path_uniprot_sprot_fasta){
  # Read gzip-compressed FASTA file:
  aa_SP <- readAAStringSet(path_uniprot_sprot_fasta)
  aa_SP <- as.data.frame(aa_SP)
  # extract the protein ID from the row-name column and append a new column with ID
  temp <- mclapply(row.names(aa_SP),
                   function(i){
                     strsplit(i, split = "|", fixed = TRUE)})
  temp <- unlist(temp)
  aa_SP$ID <- temp[seq(from = 2, to = length(temp), by=3)]
  # set meaningful headers
  colnames(aa_SP) <- c("prot_seq", "ID")
  # add column with the length of the full protein sequence
  aa_SP$prot_seq_length <- mclapply(aa_SP$prot_seq, function(i){nchar(i)})
  aa_SP$prot_seq_length <- unlist(aa_SP$prot_seq_length)
  return(aa_SP) 
}

beautifier <- function(p, x.axis.text.angle = 90, x.axis.text.hjust = NULL){
  p <- p + theme(panel.background = element_rect(fill = 'transparent', colour = NA),
                 text = element_text(),
                 legend.background = element_rect(colour = "white"),
                 legend.text = element_text(family = "sans", face='italic', hjust=0),
                 legend.key = element_rect(colour = 'white', fill = 'white'),
                 strip.background = element_rect(fill = 'transparent', colour = NA), # colour='red', fill='#CCCCFF'
                 strip.text.x = element_text(family = "sans", angle = 0),
                 strip.text.y = element_text(family = "sans", angle = 270, margin = margin(r=30)),
                 axis.text.x = element_text(family = "sans", angle = x.axis.text.angle, hjust = x.axis.text.hjust, margin=margin(1,1,2,1,"pt")),
                 axis.text.y = element_text(family = "sans", margin=margin(1,1,2,1,"pt")),
                 axis.ticks.length = unit(0.05, "cm"),
                 plot.background = element_rect(fill = 'transparent', colour = NA)
  )
  return(p)
}

paper.figure <- function(p, x.axis.text.angle = 90, x.axis.text.hjust = NULL, x.axis.text.size = 25, y.axis.text.size=25){
  p <- p + theme(panel.background = element_rect(fill = 'transparent', colour = NA),
                 text = element_text(size=25),
                 #panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                 legend.background = element_rect(colour = "white"),
                 legend.text = element_text(family = "sans", size=22, face='italic', hjust=0),
                 legend.key = element_rect(colour = 'white', fill = 'white'),
                 legend.key.size = unit(1, "cm"),
                 strip.background = element_rect(fill = 'transparent', colour = NA), # colour='red', fill='#CCCCFF'
                 strip.text.x = element_text(family = "sans",size=25, angle = 0),
                 strip.text.y = element_text(family = "sans",size=25, angle = 270, margin = margin(r=30)),
                 axis.text.x = element_text(family = "sans",size=x.axis.text.size, angle = x.axis.text.angle, hjust = x.axis.text.hjust, margin=margin(1,1,2,1,"pt")),
                 axis.text.y = element_text(family = "sans",size=y.axis.text.size, margin=margin(1,1,2,1,"pt")),
                 plot.background = element_rect(fill = 'transparent', colour = NA)
  )
  return(p)
}

tile_plot_x1x2 <- function(data, columnx1, columnx2, column, save){
  local_colors = c("white", "darkblue", "blue4", "blue3", "blue2", "blue1", "blue", "yellow","orangered3", "goldenrod2", "goldenrod1", "orangered2", "orangered", "orange", "red") #"papayawhip"
  newdata = list()
  newdata[[columnx1]] = factor(data[[columnx1]], levels = min(data[[columnx1]]):max(data[[columnx1]]))
  newdata[[columnx2]] = factor(data[[columnx2]], levels = min(data[[columnx2]]):max(data[[columnx2]]))
  newdata[[column]] = data[[column]]
  p = ggplot(as.data.frame(newdata), aes_string(x=columnx1, y=columnx2, fill=column), environment=environment())
  p = p + geom_raster() + scale_fill_gradientn(colors=parula(256), name=column, guide = "colourbar")
  p = beautifier(p)
  if (save == TRUE) {
    ggsave(paste(pathImages, column, '_x1x2', figureFormat, sep=''), width=12, height=8, dpi = 300)
  }
  return(p)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL, title="MAIN TITLE") {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout)+1, ncol(layout))))
    grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol(layout)))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row+1,
                                      layout.pos.col = matchidx$col))
    }
  }
}


## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(num_items    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$num_items)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$num_items-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

