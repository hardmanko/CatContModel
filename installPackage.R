
# Install installation helpers
install.packages(c("devtools", "R.rsp"))

# Install other dependencies from CRAN. 
install.packages(c("Rcpp", "polspline", "CircStats", "msm", "abind"))

# Install the LineChart and CMBBHT packages from github.
devtools::install_github("hardmanko/LineChart-package@v0.3.1", build_vignettes=TRUE)
devtools::install_github("hardmanko/CMBBHT@v0.1.3", build_vignettes=TRUE)

# Install this package.
#
# Installable versions of the package are attached to releases of the package.
# Their filenames are formatted as follows
#
# CatContModel_V.V.V.EXT
#
# where V.V.V is the version number and EXT is the file extension.
#
# The file extension of the file indicates what operating system it is for.
#
# Windows binary: .zip
# OSx binary:     .tgz
# Source:         .tar.gz (Any OS, but requires C++ compiler)
#
# Change PACKAGE_LOCATION to the location of the package file. 
# Change "V.V.V" to the version you want to install.
# Change the file extension (EXT) to the type of package you want to install.
#
install.packages("PACKAGE_LOCATION/CatContModel_V.V.V.EXT", repos=NULL)
