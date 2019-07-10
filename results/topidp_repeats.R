rm(list = ls(all = TRUE))
gc()
library("ggplot2")
library("plyr")
library("scales")
require(grid)
library(RColorBrewer)
library(stringr)
source("helpers.R")

sp_path = "/data/swissprot_annotations.tsv"
sp_all = load_swissprot(sp_path)

#tr_consensus_path = "results/tr_annotations/tr_consensus.csv"
#tr_consensus_path = paste(local_base_path, tr_consensus_path, sep=local_path_separator)
#tr_consensus = read.csv(tr_consensus_path, header = TRUE, sep=',')
#tr_consensus = na.omit(tr_consensus)

# tr_consensus_disorder_path = "results/tr_annotations/tr_all_with_disorder.csv"
tr_consensus_disorder_path = "results/tr_idr_overlap/repeat_disorder_overlap_regions.csv" #missing avg_idp and idr_perc
tr_consensus_disorder = load_tr_annotations(tr_consensus_disorder_path)

# Add meta_data from sp_all to tr_all. -> Do a left join.
tr_consensus_disorder_sp = merge(x = tr_consensus_disorder, y = sp_all, by = "ID", all.x = TRUE)

###################################################s############

p<-list()
pl = ggplot(tr_consensus_disorder, aes(x = avg_idp,  fill=l_type)) +
  geom_density(aes(y=..count..), alpha=0.3, size=1, position="stack") 
pl + theme_minimal() + labs(title="TOP-IDP of tandem repeats in Swissprot (stacked counts)", x="Avg Top-IDP", y="Number of repeats")

pl = ggplot(tr_consensus_disorder, aes(x = idr_perc,  fill=l_type)) +
  geom_density(aes(y=..count..), alpha=0.3, size=1, position="stack") 
pl + theme_minimal() + labs(title="% Disorder of tandem repeats in Swissprot (stacked counts)", x="Disorder % of the repeat", y="Number of repeats")

repeat_type_bin <- with(tr_consensus_disorder, cut(l_effective, breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 30, 100, 200, 500, 2000)))
tr_consensus_disorder$repeat_group = repeat_type_bin
pl <- ggplot(tr_consensus_disorder, aes(factor(repeat_group), avg_idp, fill=l_type)) 
pl + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Disorder propensity of repeats in swissprot", x="Repeat unit length", y="Average TOP-IDP of the repeat")
ggsave("figures/repeats_topidp_vs_unitlength.png")

pl <- ggplot(tr_consensus_disorder, aes(x=l_effective, y=avg_idp, colour=l_type))  + geom_point(alpha=0.5)+
  geom_smooth(method = "lm", colour="grey", se=FALSE) + xlim(0,500) +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1)) 
pl + labs(title="Disorder propensity of repeats in swissprot", x="Repeat unit length", y="Average TOP-IDP of the repeat")


pl <- ggplot(tr_consensus_disorder, aes(factor(repeat_group), idr_perc, fill=l_type)) 
pl + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  labs(title="Disorder % of repeats (Iupred Long)", x="Repeat unit length", y="Disorder %")


########################################## Split by Superkingdoms ################################################

pl = ggplot(tr_consensus_disorder_sp, aes(x = idr_perc,  colour=l_type)) +
  geom_density(size=1)  + facet_wrap(~ Superkingdom, scales = "free")
pl + theme_minimal() + labs(title="% Disorder of tandem repeats in Swissprot", x="Disorder % of the repeat", y="Density")
ggsave("figures/repeats_idrperc_density.png")


pl = ggplot(tr_consensus_disorder_sp, aes(x = avg_idp,  colour=l_type)) +
  geom_density(size=1)  + facet_wrap(~ Superkingdom, scales = "free")
pl + theme_minimal() + labs(title="TOP-IDP of tandem repeats in Swissprot", x="Avg Top-IDP", y="Density")
ggsave("figures/repeats_topidp_density.png")

repeat_type_bin <- with(tr_consensus_disorder_sp, cut(l_effective, breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 30, 100, 200, 500, 2000)))
tr_consensus_disorder_sp$repeat_group = repeat_type_bin
#pl <- ggplot(tr_consensus_disorder_sp, aes(factor(repeat_group), avg_idp, fill=l_type))  + facet_wrap(~ Superkingdom, scales = "free")
#pl + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1)) +
#  labs(title="Disorder propensity of repeats in swissprot", x="Repeat unit length", y="Average TOP-IDP of the repeat")
avg_idp_stats <- summarySE(tr_consensus_disorder_sp, measurevar="avg_idp", groupvars=c("repeat_group", "l_type", "Superkingdom"),na.rm=TRUE)
ggplot(avg_idp_stats, aes(x=repeat_group, y=avg_idp, colour=l_type)) + 
  geom_errorbar(aes(ymin=avg_idp-se, ymax=avg_idp+se), width=0.5) +
  geom_line(group = 1) +
  facet_wrap(~ Superkingdom, scales = "free") + theme_minimal() +
  labs(title="Disorder propensity of repeats in swissprot", x="Repeat unit length", y="Average TOP-IDP of the repeat") +
  geom_point()
ggsave("figures/repeats_topidp_vs_unitlength.png")

#pl <- ggplot(tr_consensus_disorder_sp, aes(factor(repeat_group), idr_perc, fill=l_type)) 
#pl + geom_boxplot() + theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1)) +  facet_wrap(~ Superkingdom, scales = "free") +
#  labs(title="Disorder % of repeats (Iupred Long)", x="Repeat unit length", y="Disorder %")
idr_perc_stats <- summarySE(tr_consensus_disorder_sp, measurevar="idr_perc", groupvars=c("repeat_group", "l_type", "Superkingdom"),na.rm=TRUE)
ggplot(idr_perc_stats, aes(x=repeat_group, y=idr_perc, colour=l_type)) + 
  geom_errorbar(aes(ymin=idr_perc-se, ymax=idr_perc+se), width=0.5) +
  geom_line(group = 1) +
  facet_wrap(~ Superkingdom, scales = "free") + theme_minimal() +
  labs(title="Disorder % of repeats (Iupred Long)", x="Repeat unit length", y="Disorder %") +
  geom_point() 

ggsave("figures/repeats_disorderperc_vs_length.png")

########################################## AA frequencies ################################################
tr_aas_path = "results/tr_annotations/tr_all_with_disorder_with_aas.csv"
tr_aas = load_tr_annotations(tr_aas_path)
repeat_type_bin <- with(tr_aas, cut(l_effective, breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 30, 100, 200, 500, 2000)))
tr_aas$repeat_group = repeat_type_bin

disorder_type_bin <- with(tr_aas, cut(idr_perc, breaks = c(-0.1,0, 0.1, 0.3, 0.5, 1.0)))
tr_aas$disorder_group = disorder_type_bin


# First disorder promoting 
disorderpromoting<-c('A','G','R', 'D', 'H', 'Q', 'S', 'K', 'E', 'P')

plot_data_column = function(data, column) {
  stats <- summarySE(data, measurevar=column, groupvars=c("repeat_group"),na.rm=TRUE)
  
  aa_freqs_path = "results/aa_freqs/aa_freqs_reformated.csv"
  aa_freqs_path = paste(local_base_path, aa_freqs_path, sep=local_path_separator)
  aa_freqs = read.csv(aa_freqs_path, header = TRUE, sep=',', row.names=1)
  
  background_value=aa_freqs[column,'background']
  p <- ggplot(stats, aes(x=stats$repeat_group, y=stats[,column])) + 
    geom_errorbar(aes(ymin=stats[,column]-stats[,'se'], ymax=stats[,column]+stats[,'se']), width=.1) +
    geom_line(group=1) +
    geom_hline(aes(yintercept=background_value), colour="#BB0000", linetype="dashed") +
    geom_point() + theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x="Repeat unit length", y=column)
}

plots <- lapply(disorderpromoting, plot_data_column, data = tr_aas)
pdf("figures/aa_freqs_in_repeats_disorderpromoting.pdf")
multiplot(plotlist=plots, cols=2, title="Frequency of disorder promoting residues")
dev.off()

# Then order promoting 
orderpromoting<-c('T','C','N','V','L','M','I','Y','F','W')
#temp <-  subset( tr_aas, select = -n ) # Remove first N column (temporary fix)
plots <- lapply(orderpromoting, plot_data_column, data = tr_aas)
# Connect all in a multiplot
pdf("figures/aa_freqs_in_repeats_orderpromoting.pdf")
multiplot(plotlist=plots, cols=2, title="Frequency of order promoting residues")
dev.off()

aln<-c()
j=1
sixdis=tr_consensus_disorder_sp[tr_consensus_disorder_sp$repeat_group=="(5,6]" & tr_consensus_disorder_sp$idr_perc>0.8 & tr_consensus_disorder_sp$Pfam_ID.y=="" ,]

for (tr in sixdis[,"MSA"]) {
  trs = strsplit(tr, " ")
  for (i in 1:length(trs)) {
    print(nchar(trs[[1]][[i]]))
    if (nchar(trs[[1]][[i]]) == 6)  {
      aln[[j]] = trs[[1]][[i]]
      j = j + 1
    }
  }
}

library("RWebLogo")
weblogo(seqs=aln, errorbars=FALSE, format='png', resolution=500, yaxis=1)


# Error bars represent standard error of the mean
library(reshape)
dfm <- melt(tr_aas[,c("repeat_group", orderpromoting, disorderpromoting)],id.vars = 1,na.rm=TRUE)
subjmeans <- cast(dfm, repeat_group~variable, mean)
df_means <- melt(subjmeans,id.vars = "repeat_group")
df_means["disorderPromoting"] <- df_means$variable %in% disorderpromoting

ggplot(df_means,aes(x = repeat_group,y = value,colour = disorderPromoting,label=variable)) +
  geom_text(aes(label=variable, size = value),hjust=0, vjust=0.5, fontface = "bold",alpha=0.9)  + scale_radius(range = c(5,25))+
  facet_wrap(~ disorderPromoting, scales = "fixed") +
  theme_minimal() +
  labs(title="Amino acid frequencies in repeats", x="Repeat unit length", y="Average AA frequency for this unit length")
#geom_bar(aes(fill = variable),position = "stack",stat="identity") 


library("RWebLogo")
aln= readLines("../fake_alignment.fasta")
weblogo(seqs=aln, errorbars=FALSE, format='png', resolution=500, size="large", color.scheme="chemistry", units="probability")


#tr_consensus_disorder_sp$OMA_ID
#tr_consensus_disorder_sp$OrthoDB_ID
pfam=read.table("../data/pfam_annotations.tab", header = TRUE, sep='\t')
tr_consensus_disorder_sp = merge(x = tr_consensus_disorder_sp, y = pfam, by = "ID", all.x = TRUE)
tr_consensus_disorder_sp$Pfam_ID <- as.character(tr_consensus_disorder_sp$Pfam_ID)

tr_consensus_disorder_sp$Pfam_ID[is.na(tr_consensus_disorder_sp$Pfam_ID)] <- as.character(tr_consensus_disorder_sp$ID[is.na(tr_consensus_disorder_sp$Pfam_ID)])

idr_perc_stats <- summarySE(tr_consensus_disorder_sp, measurevar="idr_perc", groupvars=c("repeat_group", "l_type", "Superkingdom", "Pfam_ID"),na.rm=TRUE)
#idr_perc_stats_subset <- idr_perc_stats[idr_perc_stats$num_items>10,]

idr_perc_stats_subset <-summarySE(idr_perc_stats, measurevar="idr_perc", groupvars=c("repeat_group", "l_type", "Superkingdom"),na.rm=TRUE)
ggplot(idr_perc_stats_subset, aes(x=repeat_group, y=idr_perc)) + 
  geom_errorbar(aes(ymin=idr_perc-se, ymax=idr_perc+se), width=0.5) +
  geom_line(group = 1) +
  facet_wrap(~ Superkingdom, scales = "free") + theme_minimal() +
  labs(title="Disorder % of repeats (Iupred Long)", x="Repeat unit length", y="Disorder %") +
  geom_point()


hist(idr_perc_stats[idr_perc_stats$repeat_group=="(5,6]","idr_perc"])

unitsix = idr_perc_stats[idr_perc_stats$repeat_group=="(5,6]",]
num_items_bin <- with(unitsix, cut(num_items, breaks = c(0,2, 10, 50,100,200)))
unitsix$num_items_group = num_items_bin
ggplot(unitsix, aes(idr_perc, fill=factor(num_items_group))) + geom_histogram() + 
  facet_wrap(~ Superkingdom, scales = "free") + theme_minimal() +
  labs(title="Disorder % of repeats (Iupred Long)", x="Repeat unit length", y="Disorder %") 

########################################CLustering repeats by AA profile##############################################
# Plot a SPLOM:
library(colorspace)
disorder_type_bin <- with(tr_aas, cut(idr_perc, breaks = c(-0.1,0.8, 1.0)))
tr_aas$disorder_group = disorder_type_bin

tr_aas_micro <- tr_aas[tr_aas$l_type=='micro',]
disorder_col <- rev(rainbow_hcl(2))[as.numeric(tr_aas_micro$disorder_group)]
#pairs(tr_aas_micro[,disorderpromoting], col = disorder_col,
#      lower.panel = NULL,
#      cex.labels=2, pch=19, cex = 1.2)

# Add a legend
#par(xpd = TRUE)
#legend(x = 0.05, y = 0.4, cex = 2,
#       legend = as.character(levels(tr_aas_micro$disorder_group)),
#       fill = unique(disorder_col))
#par(xpd = NA)



only_freqs = tr_aas_micro[,c(orderpromoting, disorderpromoting)]
only_freqs <- na.omit(only_freqs) # listwise deletion of missing
pc <- princomp(only_freqs)
plot(pc)
summary(pc)  
loadings(pc)
plot(pc, type='l')

# Get principal component vectors using prcomp instead of princomp
pc <- prcomp(only_freqs)
png("figures/aa_freqs_pca_disorder.png")
par(mfrow=c(2,2))
plot(pc$x[,1],pc$x[,2], col=disorder_col)
plot(pc$x[,2],pc$x[,3], col=disorder_col)
plot(pc$x[,3],pc$x[,4], col=disorder_col)
plot(pc$x[,5],pc$x[,6], col=disorder_col)
legend(x = "center",inset = 0, 
             legend = as.character(levels(tr_aas_micro$disorder_group)),
             fill = unique(disorder_col))
dev.off()

####################kmeans
wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot(only_freqs, nc=25) 

library(cluster)
#1st and 2nd component
fit <- kmeans(only_freqs, 6)
clusplot(only_freqs, fit$cluster, color=FALSE, shade=TRUE, col.p = disorder_col,
         labels=4, lines=0)

# Centroid Plot against 1st 2 discriminant functions
library(fpc)
plotcluster(only_freqs, fit$cluster)


































































































