
# Install installation helpers
install.packages(c("devtools", "R.rsp"))

# Install other dependencies from CRAN. 
install.packages(c("Rcpp", "polspline", "CircStats", "msm", "abind"))


# Install the LineChart and CMBBHT packages from github.
devtools::install_github("hardmanko/LineChart-package@v0.3.2", build_vignettes=TRUE)
devtools::install_github("hardmanko/CMBBHT@v0.1.3", build_vignettes=TRUE)


# Set the C++ compiler to use C++11 features.
# Only neccessary if you are installing from source (.tar.gz extension),
# but should be fine to run even if you are installing the .zip or .tgz versions.
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")


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
# The binary versions for Windows and OSx do not require C++ compilation.
#
# Windows binary: .zip
# OSx binary:     .tgz
# Source:         .tar.gz (Any OS, but requires C++ compiler, like RTools)
#
# Change PACKAGE_LOCATION to the location of the package file. 
# Change "V.V.V" to the version you want to install.
# Change the file extension (EXT) to the type of package you want to install.
#
install.packages("PACKAGE_LOCATION/CatContModel_V.V.V.EXT", repos=NULL)
