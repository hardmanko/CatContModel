
# Install installation helpers
install.packages(c("devtools", "R.rsp"))

# Install dependencies from CRAN. 
install.packages(c("Rcpp", "polspline", "CircStats", "msm"))

# Install the LineChart and CMBBHT packages from github.
devtools::install_github("hardmanko/LineChart-package@v0.3.0", build_vignettes=TRUE)
devtools::install_github("hardmanko/CMBBHT@v0.1.1", build_vignettes=TRUE)

# Install this package.
#
# In the github repository, compiled versions of the package 
# that are installable by R can be found in the "packaged" subfolder.
#
# If you downloaded the whole repository, change "%repoRoot%" to the root directory of the unzipped repository.
# Then change "V.V.V" to the version you want to install (typically the latest).
#
# If you downloaded just one packaged version of the package, change the path below to point to that file.
#
install.packages("%repoRoot%/packaged/CatContModel_V.V.V.zip", repos=NULL)
#
# P.S. The .zip extension indicates that it is a binary/compiled, version of the package.
# This is almost certainly what you want, not .tar.gz, which are source versions of the package
# that must be compiled with a C++ compiler.