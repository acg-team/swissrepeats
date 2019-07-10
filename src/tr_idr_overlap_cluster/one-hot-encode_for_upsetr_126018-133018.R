###############
### This file is ment to be run on a HPC cluster.
###
### Place this file ~/swissrepeat/
### create on cluster:
### /scratch/IAS/AnisGroup/delucmat/swissrepeat/results/
### /scratch/IAS/AnisGroup/delucmat/swissrepeat/data/
### /scratch/IAS/AnisGroup/delucmat/swissrepeat/figures/
### dependencies in swissrepeat directory:
### - local_config_cluster.R
### - helpers_cluster_par.R
### dependencies in data directory:
### - tr_all_overlap_par_subset.csv
### results saved in:
### /scratch/IAS/AnisGroup/delucmat/data/listinputAA_par_subset_par.csv
###############

###############
### TO DEBUG ON LOCAL MACHINE, UNCOMMENT THIS CHUNCK
### AND CONTINUE WITH ONE-HOT-ENCODING
###############
# setwd("/home/matteo/polybox/MSc_ACLS/swissrepeat/results")
# source("local_config.R")
# setwd(paste0(local_base_path,"/results"))
# 
# rm(list = ls(all = TRUE))
# gc()
# source("helpers.R")
# tr_all_overlap <- read.csv(paste0(local_base_path, "/data/tr_all_overlap.csv"), header = TRUE)
# print(paste("Tandem Repeats with disorder-overlap information are loaded"), str(tr_all_overlap))
# 
# sp_path = "/data/swissprot_annotations.tsv"
# sp_all = load_swissprot(sp_path, tr_all_overlap)
# print(paste("Swissprot Proteins are loaded", str(sp_all)))


print("###############")
print("Housekeeping")
print("###############")
# setwd("/cfs/earth/scratch/delucmat/swissrepeat/")
# source("local_config_cluster.R")
# setwd(paste0(local_base_path))

print("Set working dir to local scratch on compute node by sourcing the local config file")
source("local_config_cluster.R")
setwd(local_base_path)
print(paste("working directory set to: ", local_base_path))

rm(list = ls(all = TRUE))
gc()

source("helpers_cluster_par.R")
print("load required packages and helper functions.")

FROM <- 126018
TO <- 133018
print(paste("Subset Proteins from:", FROM, "to", TO))

print("###############")
print("Data Loading")
print("###############")
tr_all_overlap <- read.csv(paste0(local_base_path, "data/tr_all_overlap.csv"), header = TRUE)
print(paste("Tandem Repeats with disorder-overlap information are loaded"), str(tr_all_overlap))

# sp_path = "data/swissprot_annotations.tsv"
# sp_all = load_swissprot(sp_path, tr_all_overlap)
# print(paste("Swissprot Proteins are loaded", str(sp_all)))

print("###############")
print("Compute one-hot encoding for Upset plot")
print("###############")
print("summarize overlap by protein")
overlap_by_protein= ddply(tr_all_overlap, .(ID), summarize, #Aggregate disordered regions by protein, calculate disorder residues count in a protein
                          total_TR_length_prot = sum(total_repeat_length),
                          total_disorder_length_prot = sum(disordered_overlap),
                          total_overlap_prot=sum(total_overlap),
                          tail_overlap_prot=sum(tail_overlap),
                          head_overlap_prot=sum(head_overlap),
                          DisinTR_overlap_prot=sum(DisinTR_overlap),
                          TRinDis_overlap_prot=sum(TRinDis_overlap))
print(paste("Overlap regions are summarized for each TR", str(overlap_by_protein)))

# sp_overlap = merge(x = overlap_by_protein, y = sp_all, by = "ID", all.x = TRUE)
# print("meta-data from swissprot data is added to tandem repeats data", str(sp_overlap))

print("###############")
print("one-hot-encode each AA by it's set")
print("###############")
one_hot_encode_AA <- function(overlap_by_protein, FROM, TO){
  # initiallize data.frame
  listinputAA <- data.frame(matrix(ncol = 7, nrow = 0))
  
  foreach(Prot=iter(overlap_by_protein$ID[c(FROM:TO)]), .combine = 'rbind',
          # foreach(Prot=iter(overlap_by_protein$ID[which(overlap_by_protein$total_overlap_prot != 0)]), .combine = 'rbind',
          .export = c("listinputAA", "overlap_by_protein"), .verbose = TRUE) %do% {
            TR <- 0
            while (TR < overlap_by_protein$total_TR_length_prot[which(overlap_by_protein$ID == Prot)]) {
              listinputAA[(nrow(listinputAA) + 1),] <-  c(as.character(Prot),1,0,0,0,0,0)
              TR <- TR+1
            }
            Dis <- 0
            while (Dis < overlap_by_protein$total_disorder_length_prot[which(overlap_by_protein$ID == Prot)]) {
              listinputAA[(nrow(listinputAA) + 1),] <- c(as.character(Prot),0,1,0,0,0,0)
              Dis <- Dis+1
            }
            headov <- 0
            while (headov < overlap_by_protein$head_overlap_prot[which(overlap_by_protein$ID == Prot)]) {
              listinputAA[(nrow(listinputAA) + 1),] <- c(as.character(Prot),1,1,1,0,0,0)
              headov <- headov+1
            }
            tailov <- 0
            while (tailov < overlap_by_protein$tail_overlap_prot[which(overlap_by_protein$ID == Prot)]) {
              listinputAA[(nrow(listinputAA) + 1),] <- c(as.character(Prot),1,1,0,1,0,0)
              tailov <- tailov+1
            }
            DisinTRov <- 0
            while (DisinTRov < overlap_by_protein$DisinTR_overlap_prot[which(overlap_by_protein$ID == Prot)]) {
              listinputAA[(nrow(listinputAA) + 1),] <- c(as.character(Prot),1,1,0,0,1,0)
              DisinTRov <- DisinTRov+1
            }
            TRinDisov <- 0
            while (TRinDisov < overlap_by_protein$TRinDis_overlap_prot[which(overlap_by_protein$ID == Prot)]) {
              listinputAA[(nrow(listinputAA) + 1),] <- c(as.character(Prot),1,1,0,0,0,1)
              TRinDisov <- TRinDisov+1
            }
          }
  
  colnames(listinputAA) <- c("ID", "AA_in_TR", "AA_Disorder", "AA_in_headoverlap", "AA_in_tailoverlap", "AA_in_DisinTRoverlap", "AA_in_TRinDisoverlap")
  return(listinputAA)
}

start_time <- Sys.time()
print(paste("computation started at: ", start_time))

listinputAA <- one_hot_encode_AA(overlap_by_protein, FROM, TO)

end_time <- Sys.time()
run_time <- end_time-start_time
print(paste("Time used to one-hot-encode AA:", run_time))

print("###############")
print("Quick-check results")
print("###############")
print(paste("This output displays the header of the AA-encoding."))
head(listinputAA)
print(paste("This output shows the structure of the AA-encoding."))
str(listinputAA)
print(paste("Check if there are all 1"))
unique(listinputAA[-1])




print("###############")
print("Save data")
print("###############")
write.csv(listinputAA, file = paste0(local_base_path, "data/listinputAA_", FROM, "-", TO, ".csv"))
print(paste0("AA one-hot-encoding saved in: ", local_base_path, "data/listinputAA_", FROM, "-", TO, ".csv"))
