# Make a copy of the template file by entering this command in the shell:
#cp local_config_TEMPLATE.R local_config.R

# Adapt these paths to your system:
# Set the directory in which all R-packages should be installed
R_PACKAGE_PATH = "/cfs/earth/scratch/delucmat/R-packages/"

# Set the CRAN mirror of which the packages should downloaded
REPOS = "https://stat.ethz.ch/CRAN/" # more repositories: https://cran.r-project.org/mirrors.html

# local_base_path = "/cfs/earth/scratch/delucmat/swissrepeat/"
local_base_path = paste0(Sys.getenv("LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH"), "/")
print(paste0("Value obtained for local_base_path: ", local_base_path))


# *local_path_separator* is most likely '/' on a Unix system and '\' on a Windows system.
local_path_separator = "/"

# save images?
save = TRUE
pathImages = "/figures/"
figureFormat = ".png"

# get ENTREZ_Key from 
# taxize::use_entrez()
# usethis::edit_r_environ()
# and append there: ENTREZ_KEY='9892f066f3c989e29daf1c2821d9a054a108'