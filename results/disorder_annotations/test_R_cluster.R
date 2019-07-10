print(paste("Hello World from R."))

R_package_path = "/scratch/IAS/AnisGroup/delucmat/R-packages/"
.libPaths(R_package_path)

# load and if necessary install packages
load.lib<-c("plyr")
install.lib<-load.lib[!load.lib %in% installed.packages()] #check which libs are already installed
for (lib in install.lib){
  # installing/loading the libs:
  if(!require(lib)) {
    install.packages(lib, lib = R_package_path, repos = "https://stat.ethz.ch/CRAN/"); 
    require(lib, lib.loc = R_package_path) 
  } 
}
# load all the packages
lapply(load.lib, FUN = function(x){require(x, lib.loc = R_package_path)}, character.only = TRUE)
