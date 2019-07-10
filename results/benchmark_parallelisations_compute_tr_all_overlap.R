#####################################################################################
### Housekeeping
#####################################################################################
# Before executing this file, do: setwd("path/to/your/git/on/swissrepeat/results")
setwd("/home/matteo/polybox/MSc_ACLS/swissrepeat/results")
source("local_config.R")
setwd(paste0(local_base_path,"/results"))

rm(list = ls(all = TRUE))
gc()
source("helpers.R")

# Load library for benchmarking
devtools::install_github("eddelbuettel/rbenchmark")
library(rbenchmark)

#####################################################################################
### Data Loading 
#####################################################################################
tr_path = "results/tr_annotations/tr_annotations.csv"
tr_all = load_tr_annotations(tr_path)

sp_path = "data/swissprot_annotations.tsv"
sp_all = load_swissprot(sp_path, tr_all)

# Add meta_data from sp_all to tr_all. -> Do a left join.
tr_all_sp = merge(x = tr_all, y = sp_all, by = "ID", all.x = TRUE)

# Add disorder information from modiDB
discoor_path = "results/disorder_annotations/mobidb_coordinates.csv"
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

#####################################################################################
### Define different approaches w/ and w/o parallel enviroments
#####################################################################################
no_par_for <- function(tr_all, sp_all, discoor_all){
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
  
  
  for (IDprot in sp_all$ID[1:30]){ # for each protein in swissprot. Subset sp_all$ID[1:1000] takes about 10min
    if (IDprot %in% tr_all$ID){ # check if the protein has a TR
      for (IDdis in discoor_all$DIS_id[which(discoor_all$ID == IDprot)]){ # extract all disorder regions from this protein
        for (IDtr in tr_all$TR_id[which(tr_all$ID == IDprot)]){ # loop through all TRs of this protein
          if (tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
              tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
              discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) 
          {# check for each TR region if it has a tail-overlap with a disorder region
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
}

no_par_foreach_outer <- function(tr_all, sp_all, discoor_all){
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
  
  ### with foreach on the most outer loop (no paralell)
  foreach(i=1:length(sp_all$ID[1:30]), .combine = 'rbind') %do% { # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
    IDprot <- sp_all$ID[i]
    if (IDprot %in% tr_all$ID){ # check if the protein has a TR
      for (IDdis in discoor_all$DIS_id[which(discoor_all$ID == IDprot)]){ # extract all disorder regions from this protein
        for (IDtr in tr_all$TR_id[which(tr_all$ID == IDprot)]){ # loop through all TRs of this protein
          if (tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
              tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
              discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) 
          {# check for each TR region if it has a tail-overlap with a disorder region
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
}

no_par_foreach_outer_iterobj <- function(tr_all, sp_all, discoor_all){
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
  
  ### with foreach on the most outer loop and an iterator object (no paralell)
  foreach(IDprot=iter(sp_all$ID[1:30]), .combine = 'rbind') %do% { # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
    if (IDprot %in% tr_all$ID){ # check if the protein has a TR
      for (IDdis in discoor_all$DIS_id[which(discoor_all$ID == IDprot)]){ # extract all disorder regions from this protein
        for (IDtr in tr_all$TR_id[which(tr_all$ID == IDprot)]){ # loop through all TRs of this protein
          if (tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
              tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
              discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) 
          {# check for each TR region if it has a tail-overlap with a disorder region
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
}

no_par_foreach_outer_iterobj_when <- function(tr_all, sp_all, discoor_all){
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
  
  ### with foreach on the most outer loop, iterator object, when statement (no paralell)
  foreach(IDprot=iter(sp_all$ID[1:30]), .combine = 'rbind') %:%
    when(IDprot %in% tr_all$ID) %do% {   # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
      for (IDdis in discoor_all$DIS_id[which(discoor_all$ID == IDprot)]){ # extract all disorder regions from this protein
        for (IDtr in tr_all$TR_id[which(tr_all$ID == IDprot)]){ # loop through all TRs of this protein
          if (tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
              tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
              discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) 
          {# check for each TR region if it has a tail-overlap with a disorder region
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

no_par_foreach_2outer_iterobj_when <- function(tr_all, sp_all, discoor_all){
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
  
  ### with foreach on the 2 outer loops, iterator object, when statement (no paralell)
  foreach(IDprot=iter(sp_all$ID[1:30]), .combine = 'rbind') %:%
    when(IDprot %in% tr_all$ID) %do% {   # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
      foreach(IDdis=iter(discoor_all$DIS_id[which(discoor_all$ID == IDprot)]), .combine = 'rbind') %do% { # extract all disorder regions from this protein
        for (IDtr in tr_all$TR_id[which(tr_all$ID == IDprot)]){ # loop through all TRs of this protein
          if (tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
              tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
              discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) 
          {# check for each TR region if it has a tail-overlap with a disorder region
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

no_par_foreach_3outer_iterobj_when <- function(tr_all, sp_all, discoor_all){
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
  
  ### with foreach on the 3 outer loops, iterator object, when statement (no paralell)
  foreach(IDprot=iter(sp_all$ID[1:30]), .combine = 'rbind') %:%
    when(IDprot %in% tr_all$ID) %:%    # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
    foreach(IDdis=iter(discoor_all$DIS_id[which(discoor_all$ID == IDprot)]), .combine = 'rbind') %:% # extract all disorder regions from this protein
    foreach(IDtr=iter(tr_all$TR_id[which(tr_all$ID == IDprot)]), .combine = 'rbind') %do% { # loop through all TRs of this protein
      if (tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
          tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
          discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) 
      {# check for each TR region if it has a tail-overlap with a disorder region
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

par_foreach_outer <- function(tr_all, sp_all, discoor_all, numcores){
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
  
  # register cluster
  cl <- makeCluster(numcores)
  registerDoParallel(cl)
  
  ### with foreach on the most outer loop (paralell)
  foreach(i=1:length(sp_all$ID[1:30]), .combine = 'rbind') %dopar% { # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
    IDprot <- sp_all$ID[i]
    if (IDprot %in% tr_all$ID){ # check if the protein has a TR
      for (IDdis in discoor_all$DIS_id[which(discoor_all$ID == IDprot)]){ # extract all disorder regions from this protein
        for (IDtr in tr_all$TR_id[which(tr_all$ID == IDprot)]){ # loop through all TRs of this protein
          if (tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
              tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
              discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) 
          {# check for each TR region if it has a tail-overlap with a disorder region
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
  stopCluster(cl)
}

par_foreach_outer_iterobj <- function(tr_all, sp_all, discoor_all, numcores){
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
  
  # register cluster
  cl <- makeCluster(numcores)
  registerDoParallel(cl)
  
  ### with foreach on the most outer loop and an iterator object (no paralell)
  foreach(IDprot=iter(sp_all$ID[1:30]), .combine = 'rbind') %dopar% { # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
    if (IDprot %in% tr_all$ID){ # check if the protein has a TR
      for (IDdis in discoor_all$DIS_id[which(discoor_all$ID == IDprot)]){ # extract all disorder regions from this protein
        for (IDtr in tr_all$TR_id[which(tr_all$ID == IDprot)]){ # loop through all TRs of this protein
          if (tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
              tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
              discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) 
          {# check for each TR region if it has a tail-overlap with a disorder region
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
  stopCluster(cl)
}

par_foreach_outer_iterobj_when <- function(tr_all, sp_all, discoor_all, numcores){
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
  
  # register cluster
  cl <- makeCluster(numcores, type = "FORK")
  registerDoParallel(cl)
  
  ### with foreach on the first outer loop, iterator object, when statement (forking instead of parallell)
  foreach(IDprot=iter(sp_all$ID[1:300]), .combine = 'rbind', .inorder=TRUE) %:%
    when(IDprot %in% tr_all$ID) %dopar%    # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
    for (IDdis in discoor_all$DIS_id[which(discoor_all$ID == IDprot)]){ # extract all disorder regions from this protein
      for (IDtr in tr_all$TR_id[which(tr_all$ID == IDprot)]){ # loop through all TRs of this protein
        if (tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
            tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
            discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) 
        {# check for each TR region if it has a tail-overlap with a disorder region
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
  stopCluster(cl)
}

par_foreach_2outer_iterobj_when <- function(tr_all, sp_all, discoor_all, numcores){
  # # add to each TR a unique identification number: TR_id
  # tr_all <- tr_all %>% 
  #   mutate(TR_id = row_number())
  # 
  # # add to each TR the position in the protein where the TR sequence ends: end
  # tr_all <- tr_all %>% 
  #   mutate(end = begin+total_repeat_length)
  # 
  # # initallize for each TR the overlap columns
  # tr_all <- tr_all %>% 
  #   mutate(tail_overlap = 0,
  #          head_overlap = 0,
  #          DisinTR_overlap = 0,
  #          TRinDis_overlap = 0,
  #          total_overlap = 0)
  # 
  # # add to each disorder region a unique identification number: DIS_id
  # discoor_all <- discoor_all %>% 
  #   mutate(DIS_id = row_number())
  # 
  # # register cluster
  # cl <- makeCluster(numcores, type = "FORK")
  # registerDoParallel(cl)
  # 
  # ### with foreach on the 2 outer loops, iterator object, when statement (no paralell)
  # foreach(IDprot=iter(sp_all$ID[1:30]), .combine = 'rbind') %:%
  #   when(IDprot %in% tr_all$ID) %:% {   # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
  #     foreach(IDdis=iter(discoor_all$DIS_id[which(discoor_all$ID == IDprot)]), .combine = 'rbind') %dopar% { # extract all disorder regions from this protein
  #       for (IDtr in tr_all$TR_id[which(tr_all$ID == IDprot)]){ # loop through all TRs of this protein
  #         # NOTE: IDprott is not found.
  #         if (tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
  #             tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
  #             discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) 
  #         {# check for each TR region if it has a tail-overlap with a disorder region
  #           tr_all$tail_overlap[which(tr_all$TR_id == IDtr)] = tr_all$end[which(tr_all$TR_id == IDtr)] - discoor_all$start[which(discoor_all$DIS_id == IDdis)]
  #         } 
  #         if (discoor_all$end[which(discoor_all$DIS_id == IDdis)] > tr_all$begin[which(tr_all$TR_id == IDtr)] &
  #             discoor_all$end[which(discoor_all$DIS_id == IDdis)] < tr_all$end[which(tr_all$TR_id == IDtr)] &
  #             discoor_all$start[which(discoor_all$DIS_id == IDdis)] <= tr_all$begin[which(tr_all$TR_id == IDtr)])
  #         { # or if it has a head-overlap with a disorder region
  #           tr_all$head_overlap[which(tr_all$TR_id == IDtr)] = discoor_all$end[which(discoor_all$DIS_id == IDdis)] - tr_all$begin[which(tr_all$TR_id == IDtr)]
  #         } 
  #         if (tr_all$begin[which(tr_all$TR_id == IDtr)] <= discoor_all$start[which(discoor_all$DIS_id == IDdis)] & 
  #             discoor_all$end[which(discoor_all$DIS_id == IDdis)] <= tr_all$end[which(tr_all$TR_id == IDtr)])
  #         { # or if the whole disorder region lies in the TR-region
  #           tr_all$DisinTR_overlap[which(tr_all$TR_id == IDtr)] = discoor_all$disorder_region_length[which(discoor_all$DIS_id == IDdis)]
  #         } 
  #         if (tr_all$begin[which(tr_all$TR_id == IDtr)] >= discoor_all$start[which(discoor_all$DIS_id == IDdis)] & 
  #             discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)])
  #         { # or if the whole disorder region is bigger than the whole TR region. i.e. the whole TR region lies in a disorder region.
  #           tr_all$TRinDis_overlap[which(tr_all$TR_id == IDtr)] = tr_all$total_repeat_length[which(tr_all$TR_id == IDtr)] 
  #         }
  #       }
  #     }
  #   }
  # stopCluster(cl)
}

par_foreach_3outer_iterobj_when <- function(tr_all, sp_all, discoor_all, numcores){
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
  
  # register cluster
  cl <- makeCluster(numcores)
  registerDoParallel(cl)
  
  ### with foreach on the 3 outer loops, iterator object, when statement (no paralell)
  foreach(IDprot=iter(sp_all$ID[1:30]), .combine = 'rbind') %:%
    when(IDprot %in% tr_all$ID) %:%    # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
    foreach(IDdis=iter(discoor_all$DIS_id[which(discoor_all$ID == IDprot)]), .combine = 'rbind') %:% # extract all disorder regions from this protein
    foreach(IDtr=iter(tr_all$TR_id[which(tr_all$ID == IDprot)]), .combine = 'rbind') %dopar% { # loop through all TRs of this protein
      if (tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
          tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
          discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) 
      {# check for each TR region if it has a tail-overlap with a disorder region
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
  stopCluster(cl)
}


#####################################################################################
### Benchmark results
#####################################################################################
print(paste("-------------------------------------------------------------------------------------------------------"))
print("Start evaluating time for different calculation approaches.")
print("The most outer data-frame is subsetted to 30 entries. Take care when interpolating the results to the whole dataset!")
start_time <- Sys.time()
bm <- benchmark("w/o par" = no_par_for(tr_all, sp_all, discoor_all),
          "w/o par; foreach on 1. loop" = no_par_foreach_outer(tr_all, sp_all, discoor_all),
          "w/o par; foreach on 1. loop; iterator object" = no_par_foreach_outer_iterobj(tr_all, sp_all, discoor_all),
          "w/o par; foreach on 1. loop; iterator object; when statement" = no_par_foreach_outer_iterobj_when(tr_all, sp_all, discoor_all),
          # "w/o par; foreach on 1.&2. loop; iterator object" = no_par_foreach_2outer_iterobj_when(tr_all, sp_all, discoor_all),
          # "w/o par; foreach on 1.,2.&3. loop; iterator object" = no_par_foreach_3outer_iterobj_when(tr_all, sp_all, discoor_all),
          "w/ par; foreach on 1. loop" = par_foreach_outer(tr_all, sp_all, discoor_all, numcores),
          "w/ par; foreach on 1. loop; iterator object" = par_foreach_outer_iterobj(tr_all, sp_all, discoor_all, numcores),
          "w/ par; foreach on 1. loop; iterator object; when statement" = par_foreach_outer_iterobj_when(tr_all, sp_all, discoor_all, numcores),
          # "w/ par; foreach on 1.&2. loop; iterator object" = par_foreach_2outer_iterobj_when(tr_all, sp_all, discoor_all, numcores),
          # "w/ par; foreach on 1.,2.&3. loop; iterator object" = par_foreach_3outer_iterobj_when(tr_all, sp_all, discoor_all, numcores),
          replications = 1)
bm[base::order(rownames(bm), decreasing = FALSE),]
end_time <- Sys.time()
print(paste("End of benchmark. Time used for benchmarking: ", end_time-start_time))
print(paste("-------------------------------------------------------------------------------------------------------"))
# http://r.789695.n4.nabble.com/Meaning-of-proc-time-td2303263.html#a2306691
#####################################################################################
### More parallel approaches: However with difficulties.
#####################################################################################
# ### with foreach on the 3 outer loops, iterator object, multiple when statements (no paralell)
# NOTE: Multiple when statements in sequence are not evaluated correctly.
# foreach(IDprot=iter(sp_all$ID[1:30]), .combine = 'rbind') %:%
#   when(IDprot %in% tr_all$ID) %:%    # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
#   foreach(IDdis=iter(discoor_all$DIS_id[which(discoor_all$ID == IDprot)]), .combine = 'rbind') %:% # extract all disorder regions from this protein
#   foreach(IDtr=iter(tr_all$TR_id[which(tr_all$ID == IDprot)]), .combine = 'rbind') %:% # loop through all TRs of this protein
#   when(tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
#          tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
#          discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) %do%
#          {# check for each TR region if it has a tail-overlap with a disorder region
#            
#            ## This chunck is for debugging:
#            print(paste("TR-id", IDtr,
#                        # "prot-id", IDprot,
#                        # "discoor-id", IDdis,
#                        # "discoor-protid", discoor_all$DIS_id[which(discoor_all$ID == IDprot)],
#                        "TR-start", tr_all$begin[which(tr_all$TR_id == IDtr)],
#                        "TR-end", tr_all$end[which(tr_all$TR_id == IDtr)],
#                        "discorr-start", discoor_all$start[which(discoor_all$DIS_id == IDdis)],
#                        "discorr-end", discoor_all$end[which(discoor_all$DIS_id == IDdis)],
#                        "tail_overlap", tr_all$end[which(tr_all$TR_id == IDtr)] - discoor_all$start[which(discoor_all$DIS_id == IDdis)]))
#          }
# 
# ### with foreach on the 3 outer loops, iterator object, when statement (paralell)
# NOTE: Multiple when statements in sequence are not evaluated correctly.
# # register cluster
# cl <- makeCluster(numcores)
# registerDoParallel(cl)
# 
# system.time(
#   foreach(IDprot=iter(sp_all$ID[1:300]), .combine = 'rbind') %:%
#     when(IDprot %in% tr_all$ID) %:%    # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
#     foreach(IDdis=iter(discoor_all$DIS_id[which(discoor_all$ID == IDprot)]), .combine = 'rbind') %:% # extract all disorder regions from this protein
#     foreach(IDtr=iter(tr_all$TR_id[which(tr_all$ID == IDprot)]), .combine = 'rbind') %:% # loop through all TRs of this protein
#     when(tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
#            tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
#            discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) %dopar%
#            {# check for each TR region if it has a tail-overlap with a disorder region
#              
#              ## This chunck is for debugging:
#              print(paste("TR-id", IDtr,
#                          # "prot-id", IDprot,
#                          # "discoor-id", IDdis,
#                          # "discoor-protid", discoor_all$DIS_id[which(discoor_all$ID == IDprot)],
#                          "TR-start", tr_all$begin[which(tr_all$TR_id == IDtr)],
#                          "TR-end", tr_all$end[which(tr_all$TR_id == IDtr)],
#                          "discorr-start", discoor_all$start[which(discoor_all$DIS_id == IDdis)],
#                          "discorr-end", discoor_all$end[which(discoor_all$DIS_id == IDdis)],
#                          "tail_overlap", tr_all$end[which(tr_all$TR_id == IDtr)] - discoor_all$start[which(discoor_all$DIS_id == IDdis)]))
#            }
# )
# 
# stopCluster(cl)

### with foreach on the 3 outer loops, iterator object, when statement (paralell with NWS backend for chunking)
### NOTE: NWS not available for current R version. Using instead MPI
# # register cluster
# cl <- makeCluster(numcores)
# registerDoParallel(cl)
# 
# system.time(
# foreach(IDprot=iter(sp_all$ID[1:300]), .combine = 'rbind') %:%
#   when(IDprot %in% tr_all$ID) %:%    # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
#     foreach(IDdis=iter(discoor_all$DIS_id[which(discoor_all$ID == IDprot)]), .combine = 'rbind') %:% # extract all disorder regions from this protein
#       foreach(IDtr=iter(tr_all$TR_id[which(tr_all$ID == IDprot)]), .combine = 'rbind', .options.nws=list(chunkSize=2)) %:% # loop through all TRs of this protein
#         when(tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
#             tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
#             discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) %dopar%
#           {# check for each TR region if it has a tail-overlap with a disorder region
#           
#           ## This chunck is for debugging:
#           print(paste("TR-id", IDtr,
#                       # "prot-id", IDprot,
#                       # "discoor-id", IDdis,
#                       # "discoor-protid", discoor_all$DIS_id[which(discoor_all$ID == IDprot)],
#                       "TR-start", tr_all$begin[which(tr_all$TR_id == IDtr)],
#                       "TR-end", tr_all$end[which(tr_all$TR_id == IDtr)],
#                       "discorr-start", discoor_all$start[which(discoor_all$DIS_id == IDdis)],
#                       "discorr-end", discoor_all$end[which(discoor_all$DIS_id == IDdis)],
#                       "tail_overlap", tr_all$end[which(tr_all$TR_id == IDtr)] - discoor_all$start[which(discoor_all$DIS_id == IDdis)]))
#         }
# )
# 
# stopCluster(cl)


### with foreach on the 3 outer loops, iterator object, when statement (paralell MPI backend for chunking) 
### NOTE: R session crashes when closing the cluster
# # register cluster
# require(doMPI)
# cl <- startMPIcluster()
# registerDoMPI(cl)
# opt <- list(chunkSize=10)
# 
# 
# system.time(
# foreach(IDprot=iter(sp_all$ID[1:300]), .combine = 'rbind', .options.mpi=opt) %:%
#   when(IDprot %in% tr_all$ID) %:%    # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
#     foreach(IDdis=iter(discoor_all$DIS_id[which(discoor_all$ID == IDprot)]), .combine = 'rbind') %:% # extract all disorder regions from this protein
#       foreach(IDtr=iter(tr_all$TR_id[which(tr_all$ID == IDprot)]), .combine = 'rbind') %:% # loop through all TRs of this protein
#         when(tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
#             tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
#             discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) %dopar%
#           {# check for each TR region if it has a tail-overlap with a disorder region
#           
#           ## This chunck is for debugging:
#           print(paste("TR-id", IDtr,
#                       # "prot-id", IDprot,
#                       # "discoor-id", IDdis,
#                       # "discoor-protid", discoor_all$DIS_id[which(discoor_all$ID == IDprot)],
#                       "TR-start", tr_all$begin[which(tr_all$TR_id == IDtr)],
#                       "TR-end", tr_all$end[which(tr_all$TR_id == IDtr)],
#                       "discorr-start", discoor_all$start[which(discoor_all$DIS_id == IDdis)],
#                       "discorr-end", discoor_all$end[which(discoor_all$DIS_id == IDdis)],
#                       "tail_overlap", tr_all$end[which(tr_all$TR_id == IDtr)] - discoor_all$start[which(discoor_all$DIS_id == IDdis)]))
#         }
# )
# 
# closeCluster(cl)
# 
# ### TODO with foreach on the 3 outer loops, iterator object, when statement (paralell, manual chunking)
# # register cluster
# cl <- makeCluster(numcores)
# registerDoParallel(cl)
# 
# # split sp_all into equally sized chunks
# # foreach chunck of sp_all_chunk
# system.time(
#   foreach(IDprot=iter(sp_all$ID[1:300]), .combine = 'rbind') %:%
#     when(IDprot %in% tr_all$ID) %:%    # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
#     foreach(IDdis=iter(discoor_all$DIS_id[which(discoor_all$ID == IDprot)]), .combine = 'rbind') %:% # extract all disorder regions from this protein
#     foreach(IDtr=iter(tr_all$TR_id[which(tr_all$ID == IDprot)]), .combine = 'rbind') %:% # loop through all TRs of this protein
#     when(tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
#            tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
#            discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) %dopar%
#            {# check for each TR region if it has a tail-overlap with a disorder region
#              
#              ## This chunck is for debugging:
#              print(paste("TR-id", IDtr,
#                          # "prot-id", IDprot,
#                          # "discoor-id", IDdis,
#                          # "discoor-protid", discoor_all$DIS_id[which(discoor_all$ID == IDprot)],
#                          "TR-start", tr_all$begin[which(tr_all$TR_id == IDtr)],
#                          "TR-end", tr_all$end[which(tr_all$TR_id == IDtr)],
#                          "discorr-start", discoor_all$start[which(discoor_all$DIS_id == IDdis)],
#                          "discorr-end", discoor_all$end[which(discoor_all$DIS_id == IDdis)],
#                          "tail_overlap", tr_all$end[which(tr_all$TR_id == IDtr)] - discoor_all$start[which(discoor_all$DIS_id == IDdis)]))
#            }
# )
# 
# stopCluster(cl)
# 
# ### with foreach on the 3 outer loops, iterator object, when statement (forking instead of parallell)
# # register cluster
# cl <- makeCluster(numcores, type = "FORK")
# registerDoParallel(cl)
# 
# system.time(
#   foreach(IDprot=iter(sp_all$ID[1:300]), .combine = 'rbind', .inorder=TRUE) %:%
#     when(IDprot %in% tr_all$ID) %:%    # for each protein in swissprot. subset sp_all$ID[1:30] good for debugging, finds two results and takes about few seconds. Subset sp_all$ID[1:1000] takes about 10min
#     foreach(IDdis=iter(discoor_all$DIS_id[which(discoor_all$ID == IDprot)]), .combine = 'rbind', .inorder=TRUE) %:% # extract all disorder regions from this protein
#     foreach(IDtr=iter(tr_all$TR_id[which(tr_all$ID == IDprot)]), .combine = 'rbind', .inorder=TRUE) %:% # loop through all TRs of this protein
#     when(tr_all$end[which(tr_all$TR_id == IDtr)] > discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
#            tr_all$begin[which(tr_all$TR_id == IDtr)] < discoor_all$start[which(discoor_all$DIS_id == IDdis)] &
#            discoor_all$end[which(discoor_all$DIS_id == IDdis)] >= tr_all$end[which(tr_all$TR_id == IDtr)]) %dopar%
#            {# check for each TR region if it has a tail-overlap with a disorder region
#              
#              ## This chunck is for debugging:
#              print(paste("TR-id", IDtr,
#                          # "prot-id", IDprot,
#                          # "discoor-id", IDdis,
#                          # "discoor-protid", discoor_all$DIS_id[which(discoor_all$ID == IDprot)],
#                          "TR-start", tr_all$begin[which(tr_all$TR_id == IDtr)],
#                          "TR-end", tr_all$end[which(tr_all$TR_id == IDtr)],
#                          "discorr-start", discoor_all$start[which(discoor_all$DIS_id == IDdis)],
#                          "discorr-end", discoor_all$end[which(discoor_all$DIS_id == IDdis)],
#                          "tail_overlap", tr_all$end[which(tr_all$TR_id == IDtr)] - discoor_all$start[which(discoor_all$DIS_id == IDdis)]))
#            }
#     
# )
# 
# stopCluster(cl)

