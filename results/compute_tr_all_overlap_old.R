###############
### This file is ment to be run on a HPC cluster.
### Place this file ~/swissrepeat/
### create on cluster:
### /scratch/IAS/AnisGroup/delucmat/swissrepeat/results/
### /scratch/IAS/AnisGroup/delucmat/swissrepeat/data/
### /scratch/IAS/AnisGroup/delucmat/swissrepeat/figures/
### dependencies in swissrepeat directory:
### - local_config.R
### - helpers.R
### dependencies in data directory:
### - tr_annotations.csv
### - swissprot_annotations.tsv
### - modib_coordinates.csv
### results saved in:
### /scratch/IAS/AnisGroup/delucmat/results/tr_all_overlap.csv
###############


###############
### Housekeeping
###############
# setwd("/cfs/earth/scratch/delucmat/swissrepeat/")
# source("local_config_cluster.R")
# setwd(paste0(local_base_path))

# Set working dir to local scratch on compute node by sourcing the local config file
source("local_config_cluster.R")
setwd(local_base_path)

rm(list = ls(all = TRUE))
gc()
source("helpers_cluster.R")

# Data Loading 
tr_path = "data/tr_annotations.csv"
tr_all = load_tr_annotations(tr_path)

sp_path = "data/swissprot_annotations.tsv"
sp_all = load_swissprot(sp_path, tr_all)

# Add meta_data from sp_all to tr_all. -> Do a left join.
tr_all_sp = merge(x = tr_all, y = sp_all, by = "ID", all.x = TRUE)

# Add disorder information from modiDB
discoor_path = "data/mobidb_coordinates.csv"
discoor_all = load_disorder_annotations(discoor_path)

discoor_all$disorder_region_length =  (discoor_all$end - discoor_all$start) + 1
discoor_all$center = floor(discoor_all$start + (discoor_all$disorder_region_length/2))
discoor_all_sp = merge(x = discoor_all, y = sp_all, by = "ID", all.x = TRUE)

# Remove eventual regions with incorrect coordinates
discoor_all_sp = discoor_all_sp[(discoor_all_sp$center<discoor_all_sp$Length),]

#Aggregate disordered regions by protein, calculate disorder residues count in a protein
sp_protein= ddply(discoor_all, .(ID), summarize,
                  disorder_count=(sum(disorder_region_length)))
sp_protein = merge(x = sp_protein, y = sp_all, by = "ID", all.x = TRUE)

# add to each TR a unique identification number: TR_id
tr_all <- tr_all %>% 
  mutate(TR_id = row_number())

# add to each TR the position in the protein where the TR sequence ends: end
tr_all <- tr_all %>% 
  mutate(end = begin+total_repeat_length)

# initallize for each TR the overlap columns
tr_all <- tr_all %>% 
  mutate(tail_overlap = 0,
         head_overlap = 0,
         DisinTR_overlap = 0,
         TRinDis_overlap = 0,
         total_overlap = 0)

# add to each disorder region a unique identification number: DIS_id
discoor_all <- discoor_all %>% 
  mutate(DIS_id = row_number())

start_time <- Sys.time()
# determine the overlap of TR regions with disorder regions for each protein
for (IDprot in sp_all$ID[1:1000]){ # for each protein in swissprot. Subset sp_all$ID[1:1000] takes about 10min
  if (IDprot %in% tr_all$ID){ # check if the protein has a TR
    for (IDdis in discoor_all$DIS_id[which(discoor_all$ID == IDprot)]){ # extract all disorder regions from this protein
      for (IDtr in tr_all$TR_id[which(tr_all$ID == IDprot)]){ # loop through all TRs of this protein
        if (tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
            tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
            discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) 
        {# check for each TR region if it has a tail-overlap with a disorder region
          
          ### This chunck is for debugging:
          # print(paste("TR-id", IDtr,
          #             # "prot-id", IDprot,
          #             # "discoor-id", IDdis,
          #             # "discoor-protid", discoor_all$DIS_id[which(discoor_all$ID == IDprot)],
          #             "TR-start", tr_all$begin[which(tr_all$TR_id == IDtr)],
          #             "TR-end", tr_all$end[which(tr_all$TR_id == IDtr)],
          #             "discorr-start", discoor_all$start[which(discoor_all$DIS_id == IDdis)],
          #             "discorr-end", discoor_all$end[which(discoor_all$DIS_id == IDdis)],
          #             "tail_overlap", tr_all$end[which(tr_all$TR_id == IDtr)] - discoor_all$start[which(discoor_all$DIS_id == IDdis)]))
          # }}}}}
          
          tr_all$tail_overlap[which(tr_all$TR_id == IDtr)] = tr_all$end[which(tr_all$TR_id == IDtr)] - discoor_all$start[which(discoor_all$DIS_id == IDdis)]
        } 
        if (discoor_all$end[which(discoor_all$DIS_id == IDdis)] > tr_all$begin[which(tr_all$TR_id == IDtr)] &
            discoor_all$end[which(discoor_all$DIS_id == IDdis)] < tr_all$end[which(tr_all$TR_id == IDtr)] &
            discoor_all$start[which(discoor_all$DIS_id == IDdis)] <= tr_all$begin[which(tr_all$TR_id == IDtr)])
        { # or if it has a head-overlap with a disorder region
          tr_all$head_overlap[which(tr_all$TR_id == IDtr)] = discoor_all$end[which(discoor_all$DIS_id == IDdis)] - tr_all$begin[which(tr_all$TR_id == IDtr)]
        } 
        if (tr_all$begin[which(tr_all$TR_id == IDtr)] <= discoor_all$start[which(discoor_all$DIS_id == IDdis)] & 
            discoor_all$end[which(discoor_all$DIS_id == IDdis)] <= tr_all$end[which(tr_all$TR_id == IDtr)])
        { # or if the whole disorder region lies in the TR-region
          tr_all$DisinTR_overlap[which(tr_all$TR_id == IDtr)] = discoor_all$disorder_region_length[which(discoor_all$DIS_id == IDdis)]
        } 
        if (tr_all$begin[which(tr_all$TR_id == IDtr)] >= discoor_all$start[which(discoor_all$DIS_id == IDdis)] & 
            discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)])
        { # or if the whole disorder region is bigger than the whole TR region. i.e. the whole TR region lies in a disorder region.
          tr_all$TRinDis_overlap[which(tr_all$TR_id == IDtr)] = tr_all$total_repeat_length[which(tr_all$TR_id == IDtr)] 
        } 
      }
    }
  }
}
end_time <- Sys.time()
run_time <- end_time-start_time
print(paste("Time used to determine the overlap region of TRs with disorder regions:", run_time))

# check if correct
# tr_all[which(tr_all$TR_id == 172904), ] # No overlap
# discoor_all[which(discoor_all$ID == "P30443"),] # since overlap region of this protein don't fall into the TR above
# 
# tr_all[which(tr_all$TR_id == 172906), ] # example for body overlap: [1] "TR-id 172906 TR-start 335 TR-end 356 discorr-start 335 discorr-end 365 tail_overlap 21"

# summarize overlap by TR (if one TR has multiple overlap, these are aggregated here)
tr_all$total_overlap <- rowSums(tr_all[,c("tail_overlap", "head_overlap", "TRinDis_overlap", "DisinTR_overlap")])

# save tr_overlap to file:
write.csv(tr_all, file = paste0(local_base_path, "data/tr_all_overlap.csv"))
print(paste0("TR overlap region information saved in: ", local_base_path, "data/tr_all_overlap.csv"))


#### Compute one-hot encoding for Upset plot
# summarize overlap by protein
overlap_by_protein= ddply(tr_all, .(ID), summarize, #Aggregate disordered regions by protein, calculate disorder residues count in a protein
                          total_TR_length_prot = sum(total_repeat_length),
                          total_disorder_length_prot = sum(disordered_overlap),
                          total_overlap_prot=sum(total_overlap),
                          tail_overlap_prot=sum(tail_overlap),
                          head_overlap_prot=sum(head_overlap),
                          DisinTR_overlap_prot=sum(DisinTR_overlap),
                          TRinDis_overlap_prot=sum(TRinDis_overlap))
sp_overlap = merge(x = overlap_by_protein, y = sp_all, by = "ID", all.x = TRUE)

listinputAA <- data.frame(matrix(ncol = 7, nrow = 0))

# one-hot-encode each AA by it's set
start_time <- Sys.time()
for (Prot in overlap_by_protein$ID) {
  TR <- 0
  while (TR < overlap_by_protein$total_TR_length_prot[which(overlap_by_protein$ID == Prot)]) {
    listinputAA <- rbind(listinputAA, c(as.numeric(as.factor(Prot)),1,0,0,0,0,0))
    TR <- TR+1
  }
  Dis <- 0
  while (Dis < overlap_by_protein$total_disorder_length_prot[which(overlap_by_protein$ID == Prot)]) {
    listinputAA <- rbind(listinputAA, c(as.numeric(as.factor(Prot)),0,1,0,0,0,0))
    Dis <- Dis+1
  }
  headov <- 0
  while (headov < overlap_by_protein$head_overlap_prot[which(overlap_by_protein$ID == Prot)]) {
    listinputAA <- rbind(listinputAA, c(as.numeric(as.factor(Prot)),1,1,1,0,0,0))
    headov <- headov+1
  }
  tailov <- 0
  while (tailov < overlap_by_protein$tail_overlap_prot[which(overlap_by_protein$ID == Prot)]) {
    listinputAA <- rbind(listinputAA, c(as.numeric(as.factor(Prot)),1,1,0,1,0,0))
    tailov <- tailov+1
  }
  DisinTRov <- 0
  while (DisinTRov < overlap_by_protein$DisinTR_overlap_prot[which(overlap_by_protein$ID == Prot)]) {
    listinputAA <- rbind(listinputAA, c(as.numeric(as.factor(Prot)),1,1,0,0,1,0))
    DisinTRov <- DisinTRov+1
  }
  TRinDisov <- 0
  while (TRinDisov < overlap_by_protein$TRinDis_overlap_prot[which(overlap_by_protein$ID == Prot)]) {
    listinputAA <- rbind(listinputAA, c(as.numeric(as.factor(Prot)),1,1,0,0,0,1))
    TRinDisov <- TRinDisov+1
  }
}
colnames(listinputAA) <- c("ID", "AA_in_TR", "AA_Disorder", "AA_in_headoverlap", "AA_in_tailoverlap", "AA_in_DisinTRoverlap", "AA_in_TRinDisoverlap")
end_time <- Sys.time()
run_time <- end_time-start_time
print(paste("Time used to one-hot-encode AA:", run_time))

# save listinputAA to file:
write.csv(listinputAA, file = paste0(local_base_path, "data/listinputAA.csv"))
print(paste0("AA one-hot-encoding saved in: ", local_base_path, "data/listinputAA.csv"))
