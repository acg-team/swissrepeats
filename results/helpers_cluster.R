###############
### All helper-functions for the cluster.
### Compared to the local helper file, here are only the packages loaded, which are really needed.
### If they are not yet installed, they are going to be installed.
### and the local_conig_cluster.R is loaded instead of local_config.R
###############

# require("plyr")
# require("pals") # nice color gradients
# require("ggplot2")
# require("scales")
# require("grid")
# require("RColorBrewer")
# require("tidyverse")
# require("ggrepel")
# require("xtable")
# require("Biostrings")
# require("doParallel")
# require("foreach")
# require("taxize")



# In the local config file, *local_base_path* is defined. *local_base_path* is the base path where 
# all results files that need to be read in are located.
# In the local config file, *local_path_separator* is defined. *local_path_separator* is most likely
# '/' on a Unix system and '\' on a Windows system.
source("local_config_cluster.R")

### Installing and loading R packages on Cluster
# Set all the requested libraries
request.lib<-c("plyr", "tidyverse")

# # Create a directory where R packages are installed
# if (!dir.exists(R_PACKAGE_PATH)){
#   dir.create(R_PACKAGE_PATH)
# } else {
#   print(paste(R_PACKAGE_PATH, " already exists!"))
# }
# 
# # Add the personal R-package directory to the list of which R is searching for libraries
# .libPaths(R_PACKAGE_PATH)


### Installation and Loading of libraries
# install and then load the requested packages if not already installed
install.lib<-request.lib[!request.lib %in% installed.packages()] #check which libs are already installed
for (lib in install.lib){
  # installing/loading the libs:
  if(!require(lib)) {
    install.packages(lib, repos = REPOS);
    # require(lib)
  }
}

# load requested packages if already installed
load.lib <- request.lib[!request.lib %in% install.lib]
# lapply(load.lib, require, character.only = TRUE)
lapply(request.lib, require, character.only = TRUE)

#### Task specific functions
# Order taken from Uversky's paper: http://www.tandfonline.com/doi/full/10.4161/idp.24684
# Not mentioned in Uversky's paper: "B", "O", "U", "Z", "X". These guys might need to fit in with the rest (if possible, as some of them represent multiple aa.)
aa_order_promoting_to_disorder_promoting = c("C", "W", "I", "Y", "F", "L", "H", "V", "N", "M", "R", "T", "D", "G", "A", "K", "Q", "S", "E", "P", "B", "O", "U", "Z", "X")

# Not in operator
'%!in%' <- function(x,y)!('%in%'(x,y))

# detect available cores (DO NOT RUN ON CLUSTER)
#numcores <- detectCores()-1


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


beautifier <- function(p){
  p <- p + theme(panel.background = element_rect(fill = 'transparent', colour = NA),
                 text = element_text(),
                 legend.background = element_rect(colour = "white"),
                 legend.text = element_text(family = "sans", face='italic', hjust=0),
                 legend.key = element_rect(colour = 'white', fill = 'white'),
                 strip.background = element_rect(fill = 'transparent', colour = NA), # colour='red', fill='#CCCCFF'
                 strip.text.x = element_text(family = "sans", angle = 0),
                 strip.text.y = element_text(family = "sans", angle = 270, margin = margin(r=30)),
                 axis.text.x = element_text(family = "sans", angle = 90, margin=margin(1,1,2,1,"pt")),
                 axis.text.y = element_text(family = "sans", margin=margin(1,1,2,1,"pt")),
                 axis.ticks.length = unit(0.05, "cm"),
                 plot.background = element_rect(fill = 'transparent', colour = NA)
  )
  return(p)
}

paper.figure <- function(p){
  p <- p + theme(panel.background = element_rect(fill = 'transparent', colour = NA),
                 text = element_text(size=25),
                 #panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                 legend.background = element_rect(colour = "white"),
                 legend.text = element_text(family = "sans", size=15, face='italic', hjust=0),
                 legend.key = element_rect(colour = 'white', fill = 'white'),
                 strip.background = element_rect(fill = 'transparent', colour = NA), # colour='red', fill='#CCCCFF'
                 strip.text.x = element_text(family = "sans",size=15, angle = 0),
                 strip.text.y = element_text(family = "sans",size=15, angle = 270, margin = margin(r=30)),
                 axis.text.x = element_text(family = "sans",size=15, angle = 90, margin=margin(1,1,2,1,"pt")),
                 axis.text.y = element_text(family = "sans",size=15, margin=margin(1,1,2,1,"pt")),
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

