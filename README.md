# swissrepeat

Collaboration with ZHAW and Stockholm University on annotation and analysis of the SwissProt database in terms of tandem repeats and disorder.

## Publication

Delucchi, M.; Schaper, E.; Sachenkova, O.; Elofsson, A.; Anisimova, M.  
**A New Census of Protein Tandem Repeats and Their Relationship with Intrinsic Disorder.**  
Genes 2020, 11, 407. <https://doi.org/10.3390/genes11040407>

# Getting Started

These instructions will get you a copy of the project up and running on your local machine to reproduce the published results.

## Prerequisites

```
R version 3.6.2
R studio version >1.2.5
```

## Installing
Clone this repository and cd in ```./results/```  

Make a copy of ```local_config_TEMPLATE.R``` by:
```
cp local_config_TEMPLATE.R local_config.R
```

In the newly created ```local_config.R``` specify system specific configurations.

Install the required R-packages:
```
install.packages("plyr")
install.packages("pals") # nice color gradients
install.packages("ggplot2")
install.packages("scales")
install.packages("grid")
install.packages("RColorBrewer")
install.packages("tidyverse")
install.packages("ggrepel")
install.packages("xtable")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
install.packages("doParallel")
install.packages("foreach")
install.packages("taxize")
install.packages("UpSetR")
install.packages("reshape2")
```

