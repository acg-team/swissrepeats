# Make a copy of the template file by entering this command in the shell:
#cp local_config_TEMPLATE.R local_config.R

# Adapt these paths to your system:
local_base_path = "/home/matteo/polybox/MSc_ACLS/swissrepeat"


# *local_path_separator* is most likely '/' on a Unix system and '\' on a Windows system.
local_path_separator = "/"

# save images?
save = TRUE
pathImages = "/home/matteo/polybox/MSc_ACLS/swissrepeat/paper/figures/"
figureFormat = ".png"

# get ENTREZ_Key from 
# taxize::use_entrez()
# usethis::edit_r_environ()
# and append there: ENTREZ_KEY='9892f066f3c989e29daf1c2821d9a054a108'