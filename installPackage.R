
# Install installation helpers
install.packages(c("devtools", "R.rsp"))

# Install dependencies from CRAN. 
install.packages(c("Rcpp", "polspline", "CircStats", "msm"))

# Install the LineChart and CMBBHT packages from github.
devtools::install_github("hardmanko/LineChart-package@v0.3.0", build_vignettes=TRUE)
devtools::install_github("hardmanko/CMBBHT@v0.1.1", build_vignettes=TRUE)

# Install this package.
#
# Installable versions of the package can be found in the "packaged" subfolder
# of the git repository.
#
# Depending on what operating system you have, you should use different files.
#
# On Windows, you can install precompiled versions of the package that end in ".zip".
# The precompiled package will only work on Windows.
#
# On Linux or OSx, you can install source versions of the package that end in ".tar.gz".
# You may be able to install source packages on Windows, but you need extra stuff (RTools).
#
# Change PACKAGE_LOCATION to the location of the package file. 
# Change "V.V.V" to the version you want to install (typically the latest).
# Change the extension to the type of package you want to install.

install.packages("PACKAGE_LOCATION/CatContModel_V.V.V.zip", repos=NULL)
