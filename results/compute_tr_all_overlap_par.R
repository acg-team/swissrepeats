###############
### This file is ment to be run on a HPC cluster.
### IT IS NOT RUNNING WITH A PARALLEL FRAMEWORK.
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
### - tr_annotations.csv
### - swissprot_annotations.tsv
### - modib_coordinates.csv
### results saved in:
### /scratch/IAS/AnisGroup/delucmat/results/tr_all_overlap.csv
###############

###############
### TO DEBUG ON LOCAL MACHINE, UNCOMMENT THIS CHUNCK
### AND CONTINUE WITH PREPROCESSING
###############
# setwd("/home/matteo/polybox/MSc_ACLS/swissrepeat/results")
# source("local_config.R")
# setwd(paste0(local_base_path,"/results"))
# 
# rm(list = ls(all = TRUE))
# gc()
# source("helpers.R")
# tr_path = "results/tr_annotations/tr_annotations.csv"
# sp_path = "data/swissprot_annotations.tsv"
# discoor_path = "results/disorder_annotations/mobidb_coordinates.csv"
# tr_all = load_tr_annotations(tr_path)
# sp_all = load_swissprot(sp_path, tr_all)
# discoor_all = load_disorder_annotations(discoor_path)


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

print("###############")
print("Data Loading")
print("###############")
tr_path = "data/tr_annotations.csv"
tr_all = load_tr_annotations(tr_path)
print(paste("Tandem Repeats are loaded"), str(tr_all))

sp_path = "data/swissprot_annotations.tsv"
sp_all = load_swissprot(sp_path, tr_all)
print(paste("Swissprot Proteins are loaded", str(sp_all)))

discoor_path = "data/mobidb_coordinates.csv"
discoor_all = load_disorder_annotations(discoor_path)
print(paste("Disorder coordinates are loaded", str(discoor_all)))

print("###############")
print("Data Preprocessing")
print("###############")
tr_all_sp = merge(x = tr_all, y = sp_all, by = "ID", all.x = TRUE)
print(paste("meta-data from swissprot data is added to tandem repeats data", str(tr_all_sp)))

print("Add disorder information from modiDB")
discoor_all$disorder_region_length =  (discoor_all$end - discoor_all$start) + 1
discoor_all$center = floor(discoor_all$start + (discoor_all$disorder_region_length/2))
discoor_all_sp = merge(x = discoor_all, y = sp_all, by = "ID", all.x = TRUE)

print("Remove eventual regions with incorrect coordinates")
discoor_all_sp = discoor_all_sp[(discoor_all_sp$center<discoor_all_sp$Length),]

# print("Aggregate disordered regions by protein, calculate disorder residues count in a protein")
# sp_protein= ddply(discoor_all, .(ID), summarize,
#                   disorder_count=(sum(disorder_region_length)))
# sp_protein = merge(x = sp_protein, y = sp_all, by = "ID", all.x = TRUE)
# 
# print("add to each TR a unique identification number: TR_id")
# tr_all <- tr_all %>% 
#   mutate(TR_id = row_number())

print("add to each TR the position in the protein where the TR sequence ends: end")
tr_all <- tr_all %>% 
  mutate(end = begin+total_repeat_length)
# 
# print("initallize for each TR an overlap column")
# tr_all <- tr_all %>% 
#   mutate(tail_overlap = 0,
#          head_overlap = 0,
#          DisinTR_overlap = 0,
#          TRinDis_overlap = 0,
#          total_overlap = 0)

print("add to each disorder region a unique identification number: DIS_id")
discoor_all <- discoor_all %>% 
  mutate(DIS_id = row_number())

print("###############")
print("Find disorder-overlaps in tandem repeats")
print("###############")
#### new version (should be faster):
# tr_all[which(tr_all$ID == "B4DXR9"),-2]
# discoor_all[which(discoor_all$ID == "B4DXR9"),]
# sp_all[which(sp_all$ID == "B4DXR9"),]

require(data.table)
# start and end should be labeled the same in tr_all and discoord all
names(tr_all)[names(tr_all) == "begin"] <- "start"

# foverlap works only with data.tables
tr_all <- as.data.table(tr_all)
sp_all <- as.data.table(sp_all) 
discoor_all <- as.data.table(discoor_all)

# subset the columns which are relevant for the LUT: ID, start, end
tr_all2 <- tr_all[,c(1,3,15)]
discoor_all2 <- discoor_all[,c(1,2,3)]

setkey(discoor_all2, ID, start, end)
setkey(tr_all2, ID, start, end)

# initiallize dataframe to store the results
tr_all_overlap <- tr_all

#### DisinTR
over <- foverlaps(y = tr_all2, 
          x = discoor_all2, 
          type = c("within"))
# check
# over[which(over$ID == "B4DXR9"),]
DisinTR <- over[which(!is.na(over$start) & !is.na(over$end)),] %>% # take only overlapping regions
  group_by(ID) %>% 
  mutate(diff = c(i.end - i.start)) %>% # extract the overlap range for each overlap
  mutate(DisinTR_overlap = cumsum(diff)) %>% # if there are multiple overlaps of the same kind,
  filter(DisinTR_overlap == max(DisinTR_overlap)) # add them up and keep the final value.
# check
# DisinTR[which(DisinTR$ID == "B4DXR9"),]
DisinTR <- DisinTR[,c(1,7)]
# full left join to tr_all_overlap
tr_all_overlap <- merge(x = tr_all_overlap, y = DisinTR, by = "ID", all = TRUE)
# replace the NAs with 0s
tr_all_overlap$DisinTR_overlap[is.na(tr_all_overlap$DisinTR_overlap)] <- 0
# check
# tr_all_overlap[which(tr_all_overlap$ID == "B4DXR9"), -2]



#### TRinDis
over <- foverlaps(x = tr_all2, 
                  y = discoor_all2, 
                  type = c("within"))
# check
# over[which(over$ID == "A0A023PXK7"),]
TRinDis <- over[which(!is.na(over$start) & !is.na(over$end)),] %>% # take only overlapping regions
  group_by(ID) %>% 
  mutate(diff = c(i.end - i.start)) %>% # extract the overlap range for each overlap
  mutate(TRinDis_overlap = cumsum(diff)) %>% # if there are multiple overlaps of the same kind,
  filter(TRinDis_overlap == max(TRinDis_overlap)) # add them up and keep the final value.
# check
# TRinDis[which(TRinDis$ID == "A0A023PXK7"),]
TRinDis <- as.data.frame(TRinDis[,c(1,7)])
# full left join to tr_all_overlap
tr_all_overlap <- merge(x = tr_all_overlap, y = TRinDis, by = "ID", all = TRUE)
# replace the NAs with 0s
tr_all_overlap$TRinDis_overlap[is.na(tr_all_overlap$TRinDis_overlap)] <- 0
# check
# tr_all_overlap[which(tr_all_overlap$ID == "A0A023PXK7"), -2]
# tr_all[which(tr_all$ID == "A0A023PXK7"),-2]
# discoor_all[which(discoor_all$ID == "A0A023PXK7"),]




#### head_overlap
# NOTE: head and tail overlap has a different definition in foverlaps than required. Hence it doesn't work.
# tr_all_overlap_old <- read.csv(paste0(local_base_path, local_path_separator, "data", local_path_separator, "tr_all_overlap.csv"), header = TRUE)
over <- foverlaps(y = tr_all2, 
                  x = discoor_all2, 
                  type = c("start"))

# ## experiment:
# x = data.table(chr=c("Chr1", "Chr1", "Chr2", "Chr2", "Chr2"),
#                start=c(5,10, 1, 25, 50), end=c(11,20,4,52,60))
# y = data.table(chr=c("Chr1", "Chr1", "Chr2"), start=c(1, 15,1),
#                end=c(4, 18, 55), geneid=letters[1:3])
# setkey(y, chr, start, end)
# foverlaps(x, y, type="start")
# ### end of experiment
# check
over[which(over$ID == "A0A023PXQ4"),]
over[which(over$ID == "A0A357"),]

headover <- over[which(!is.na(over$start) & !is.na(over$end)),] %>% # take only overlapping regions
  group_by(ID) %>% 
  mutate(diff = c(i.end - i.start)) %>% # extract the overlap range for each overlap
  mutate(head_overlap = cumsum(diff)) %>% # if there are multiple overlaps of the same kind,
  filter(head_overlap == max(head_overlap)) # add them up and keep the final value.
# check
headover[which(headover$ID == "A0A023PXQ4"),]
headover <- as.data.frame(headover[,c(1,7)])
# full left join to tr_all_overlap
tr_all_overlap <- merge(x = tr_all_overlap, y = headover, by = "ID", all = TRUE)
# replace the NAs with 0s
tr_all_overlap$head_overlap[is.na(tr_all_overlap$head_overlap)] <- 0
# check
# TODO: why is the 3 disorder-overlap not detected???
tr_all_overlap[which(tr_all_overlap$ID == "Q20CR3"), -2]
tr_all[which(tr_all$ID == "Q20CR3"),-2]
discoor_all[which(discoor_all$ID == "Q20CR3"),]
tr_all_overlap[which(tr_all_overlap$ID == "A0A023PXQ4"), -2]
tr_all[which(tr_all$ID == "A0A023PXQ4"),-2]
discoor_all[which(discoor_all$ID == "A0A023PXQ4"),]



#### tail_overlap
# tr_all_overlap_old <- read.csv(paste0(local_base_path, local_path_separator, "data", local_path_separator, "tr_all_overlap.csv"), header = TRUE)
over <- foverlaps(x = tr_all2, 
                  y = discoor_all2, 
                  type = c("end"))
# check
over[which(over$ID == "A0A357"),]
tailover <- over[which(!is.na(over$start) & !is.na(over$end)),] %>% # take only overlapping regions
  group_by(ID) %>% 
  mutate(diff = c(i.end - i.start)) %>% # extract the overlap range for each overlap
  mutate(tail_overlap = cumsum(diff)) %>% # if there are multiple overlaps of the same kind,
  filter(tail_overlap == max(tail_overlap)) # add them up and keep the final value.
# check
tailover[which(tailover$ID == "A0A357"),]
tailover <- as.data.frame(tailover[,c(1,7)])
# full left join to tr_all_overlap
tr_all_overlap <- merge(x = tr_all_overlap, y = tailover, by = "ID", all = TRUE)
# replace the NAs with 0s
tr_all_overlap$tail_overlap[is.na(tr_all_overlap$tail_overlap)] <- 0
# check
# TODO: why is the 3 disorder-overlap not detected???
tr_all_overlap[which(tr_all_overlap$ID == "A0A357"), -2]
tr_all[which(tr_all$ID == "A0A357"),-2]
discoor_all[which(discoor_all$ID == "A0A357"),]

#### TODO
## fix the head and tail overlap
## rerun ven diagramm
## rerun one-hot-encoding



#### Previous version:
start_time <- Sys.time()
print(paste("computation started at: ", start_time))

# determine the overlap of TR regions with disorder regions for each protein
foreach(IDprot=iter(sp_all$ID[1:30]), .combine = 'rbind', .inorder=TRUE, # BUG: using parallel backend outputs NULL. 
        .export = c("sp_all", "tr_all", "discoor_all"),
        .verbose=FALSE) %do% { 
          # for each protein in swissprot. 
  if (IDprot %in% tr_all$ID){ # check if the protein has a TR
    for (IDdis in discoor_all$DIS_id[which(discoor_all$ID == IDprot)]){ # extract all disorder regions from this protein
      for (IDtr in tr_all$TR_id[which(tr_all$ID == IDprot)]){ # loop through all TRs of this protein
        if (tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
            tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
            discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) 
        {# check for each TR region if it has a tail-overlap with a disorder region
          
          # ## This chunck is for debugging:
          # ## In the first for-loop, subset sp_all$ID[1:30] which works good for debugging -> finds two results and takes about few seconds. 
          # ## Subset sp_all$ID[1:1000] takes about 10min.
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

print("###############")
print("post-processing")
print("###############")
print("summarize overlap by TR (if one TR has multiple overlap, these are aggregated here")
tr_all$total_overlap <- rowSums(tr_all[,c("tail_overlap", "head_overlap", "TRinDis_overlap", "DisinTR_overlap")])

print("###############")
print("Quick-check results")
print("###############")
# ## FOR DEBUGGING, UNCOMMENT. Load file from calculated on cluster:
# tr_all_overlap <- read.csv(paste0(local_base_path, "/data/tr_all_overlap_par_subset.csv"), header = TRUE)

print("In this Protein which has a TR should be no overlap of the TR-region with the discoor-region...")
tr_all[which(tr_all$TR_id == 172904), ] # No overlap

print("... because the discoor-region doesn't fall into its TR-region.")
discoor_all[which(discoor_all$ID == "P30443"),] # since overlap region of this protein don't fall into the TR above

print("Here should be a body-overlap (TRinDis_overlap) of 21 AA beeing detected:")
tr_all[which(tr_all$TR_id == 172906), ] # example for body overlap: [1] "TR-id 172906 TR-start 335 TR-end 356 discorr-start 335 discorr-end 365 tail_overlap 21"

print("###############")
print("Save data")
print("###############")
write.csv(tr_all, file = paste0(local_base_path, "data/tr_all_overlap_par_subset.csv"))
print(paste0("TR overlap region information saved in: ", local_base_path, "data/tr_all_overlap_par_subset.csv"))

print("###############")
print("Compute one-hot encoding for Upset plot")
print("###############")
print("summarize overlap by protein")
overlap_by_protein= ddply(tr_all, .(ID), summarize, #Aggregate disordered regions by protein, calculate disorder residues count in a protein
                          total_TR_length_prot = sum(total_repeat_length),
                          total_disorder_length_prot = sum(disordered_overlap),
                          total_overlap_prot=sum(total_overlap),
                          tail_overlap_prot=sum(tail_overlap),
                          head_overlap_prot=sum(head_overlap),
                          DisinTR_overlap_prot=sum(DisinTR_overlap),
                          TRinDis_overlap_prot=sum(TRinDis_overlap))
print(paste("Overlap regions are summarized for each TR", str(overlap_by_protein)))

sp_overlap = merge(x = overlap_by_protein, y = sp_all, by = "ID", all.x = TRUE)
print("meta-data from swissprot data is added to tandem repeats data", str(sp_overlap))

print("###############")
print("one-hot-encode each AA by it's set")
print("###############")
# initiallize data.frame
listinputAA <- data.frame(matrix(ncol = 7, nrow = 0))

start_time <- Sys.time()
print(paste("computation started at: ", start_time))

foreach(Prot=iter(overlap_by_protein$ID[1:30]), .combine = 'rbind',
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
write.csv(listinputAA, file = paste0(local_base_path, "data/listinputAA_par_subset_par.csv"))
print(paste0("AA one-hot-encoding saved in: ", local_base_path, "data/listinputAA_par_subset_par.csv"))
