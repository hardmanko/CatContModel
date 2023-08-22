
# This file is useful to developers of the CatContModel package, not so much for users of the package.
# Users should see "installPackage.R" to install the package.

library(Rcpp)
library(devtools)
library(roxygen2)
library(usethis)

# Configuration options
packageName = "CatContModel"
CatContPackageVersion = "0.9.0"
addingDataSets = FALSE


baseDir = "D:/Programming/R/CatContModel/"

# Location of the properly formatted R package version of CatContModel
packageLocation = paste0(baseDir, packageName)

# Things now force the assumption that you are working in the package dir. Progress!
setwd(packageLocation)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # Compile using C++11
Sys.setenv("MAKEFLAGS"="-j4") # Multithreaded compilation of C++ files




######
# Select files to be copied into the package

modelFilesDir = paste0(baseDir, "source/cpp/")

cppHS = function(name) {
  c(paste0(name, ".h"), paste0(name, ".cpp"))
}


# Basic files (not in subdir)
generalCppFiles = c("Compilation.h", "CCM_Main.h", 
               cppHS("R_Interface"),
               cppHS("CCM_Data"),
               cppHS("CCM_ModelConfig"),
               cppHS("CCM_Util"), 
               cppHS("CCM_Weights"), 
               cppHS("CCM_Circular"), 
               cppHS("CCM_Linear"),
               cppHS("CCM_EqualityConstraints"), 
               cppHS("CCM_DistributionLUTs")
               )


# MemFirst files
memFirstFiles = c(cppHS("MF_RunModel"), "MF_RunModel_Misc.cpp",
                  cppHS("MF_Parameters"),
                  cppHS("MF_ModelUtil"),
                  cppHS("MF_Likelihood"))

memFirstFiles = paste0("MemFirst/", memFirstFiles)


# CatFirst files
catFirstFiles = c(cppHS("CF_Parameters"),
                  cppHS("CF_Likelihood"),
                  cppHS("CF_ModelConfig"),
                  cppHS("CF_Model"))

catFirstFiles = paste0("CatFirst/", catFirstFiles)



# Gibbs sampler files
gsFiles = c(cppHS("GibbsSampler"), 
            cppHS("GibbsParameters"), 
            cppHS("UtilityFunctions")
            )

gsFiles = paste0("GibbsSampler/", gsFiles)


cpp_files_brief = c(generalCppFiles, memFirstFiles, catFirstFiles, gsFiles)

cpp_files = paste0(modelFilesDir, cpp_files_brief)


# Get R filenames
rFilesDir = paste0(baseDir, "/source/R/")
rFiles = dir(rFilesDir)
rFiles = paste0(rFilesDir, rFiles)



# Include the package documentation
rFiles = c(rFiles, paste0(baseDir, "docs/toCopy/CatContModel-package.R") )
# And the data set documentation
if (addingDataSets) {
  rFiles = c(rFiles, paste0(baseDir, "examples/betweenItem/betweenDataDoc.R") )
}




# Delete the old files
packageFiles = dir( packageLocation, recursive = TRUE)
file.remove( packageFiles )

# Move new files
Rcpp.package.skeleton(packageName, path=baseDir, 
                      cpp_files=cpp_files, code_files=rFiles, 
                      force=TRUE, example_code=FALSE, 
                      author="Kyle O Hardman", 
                      email="kylehardman@gmail.com",
                      license="MIT + file LICENSE")


# Data sets
if (addingDataSets) {
  betweenItemData = read.delim( paste0(baseDir, "examples/betweenItem/betweenItemData.txt") )
  devtools::use_data(betweenItemData, pkg = packageLocation)
  betweenItemParameters = read.delim( paste0(baseDir, "examples/betweenItem/betweenItemParameters.txt") )
  devtools::use_data(betweenItemParameters, pkg = packageLocation)
}


###########################
# Append Imports and Suggests to DESCRIPTION

usethis::use_package("CircStats")
usethis::use_package("polspline")
usethis::use_package("msm")
usethis::use_package("abind")
usethis::use_package("stringr")

# My packages from github
usethis::use_package("LineChart")
usethis::use_package("CMBBHT")

#usethis::use_dev_package("LineChart", remote="hardmanko/LineChart-package@v0.3.2")
#usethis::use_dev_package("CMBBHT", remote="hardmanko/CMBBHT@v0.1.4")

usethis::use_package("parallel", type="Suggests")
usethis::use_package("R.rsp", type = "Suggests")
# OR use Introduction.Rnw


###########################
# Modify the DESCRIPTION

addDcfElement = function(dcf, name, value) {
  m = matrix(value, nrow=1, ncol=1, dimnames = list(c(), name))
  cbind(dcf, m)
}

dcf = read.dcf("DESCRIPTION")
dcf[,"Title"] = "Categorical and Continuous Working Memory Models for Delayed-Estimation Tasks"
dcf[,"Description"] = "Perform parameter estimation, posterior distribution analysis, and model comparison with the models used by Hardman, Vergauwe, & Ricker (2017). The models in this package are for delayed-estimation tasks that are commonly used in the working memory literature. The models are difficult to implement and work with for a variety of reasons, hence the value of this package. Hierarchical Bayesian implementations of between-item and within-item model variants used by Hardman, Vergauwe, and Ricker are included, as is the Zhang & Luck (2008) model. For any of these models, functions in this package allow you to relatively easily estimate the model parameters, plot parameter values, calculate posterior means and credible intervals, perform tests of the effect of task conditions on parameters, and calculate model fit statistics, among other things."
dcf[,"Author"] = "Kyle O. Hardman"
dcf[,"Maintainer"] = "Kyle O. Hardman <kylehardman@gmail.com>"
dcf[,"Version"] = CatContPackageVersion


dcf = addDcfElement(dcf, "Roxygen", "list(markdown = TRUE)")
dcf = addDcfElement(dcf, "Depends", "R (>= 4.1)")
dcf = addDcfElement(dcf, "VignetteBuilder", "R.rsp")

write.dcf(dcf, file="DESCRIPTION")


##############
# Get rid of the autogenerated package documentation
autogenToRM = c("man/CatContModel-package.Rd", "Read-and-delete-me") #"NAMESPACE"
file.remove( autogenToRM )


# Try these if package claims to be in use:
# unloadNamespace(packageName)
# remove.packages("CatContModel")
# library(CatContModel) # Then restart R again


######
# Basic install without documentation and with all functions exported
devtools::install(args="--no-multiarch", dependencies = FALSE)



# For basic build and install, stop here
########################################




#####################
# Quick and dirty update files and recompile
# Doesn't properly re-install the package if it is already loaded

if (FALSE) {
  
  quickAndDirtyRebuild = function() {
    
    toFilesSplit = c()
    for (mf in cpp_files_brief) {
      if (grepl("/", mf, fixed=TRUE)) {
        mf = strsplit(mf, "/", fixed=TRUE)[[1]][2]
      }
      toFilesSplit = c(toFilesSplit, mf)
    }
    
    file.copy(from=paste0(modelFilesDir, cpp_files_brief), 
              to=paste0("src/", toFilesSplit),
              overwrite=TRUE
    )
    
    if ("package:CatContModel" %in% search()) {
      detach("package:CatContModel", character.only = TRUE)
    }
    
    devtools::install(args="--no-multiarch", dependencies = FALSE)
    
  }
  
  ###
  quickAndDirtyRebuild()
  ###
}




######
# Install with documentation which also hides unexported functions

devtools::document() # If you remove NAMESPACE before documenting, the package is invalid.
file.remove("NAMESPACE") # roxygen2 no longer overwrites NAMESPACE, so remove it.
devtools::document() # If NAMESPACE is removed after 1 documentation pass, the package is valid.

# Install with documentation
devtools::install(args="--no-multiarch")




###################################################################
# Copy over LICENSE
dir.create("inst/")
dir.create("inst/doc/")
file.copy(from = paste0(baseDir, "LICENSE.md"), to="LICENSE")

# Copy over manual and supporting asis file
file.copy(from = paste0(baseDir, "docs/introduction/Introduction.pdf"), 
          to = "inst/doc/Introduction.pdf")

file.copy(from = paste0(baseDir, "docs/toCopy/Introduction.pdf.asis"), 
          to = c("inst/doc/Introduction.pdf.asis"  #, "vignettes/Introduction.pdf.asis"
          )
)



#####
# At this point, there should be a proper R package structure in packageLocation
#####

unloadNamespace(packageName)

# Remove object files cuz it's confuzzling to gcc (how does gcc not know what an object file is?)
srcDir = paste0(packageLocation, "/src/")
objectFiles = dir( srcDir, recursive = FALSE)
objectFiles = objectFiles[ grepl("\\.o$", objectFiles) ]
file.remove( paste0(srcDir, objectFiles) )

# Install the package from source. This is the full installation.
install.packages(pkgs=packageLocation, repos=NULL, type="source", 
                 build_vignettes=TRUE, clean=FALSE)



#################################################
# Use the above for simple test installations.  #
# Use the below for package publication.        #
#                                               #
# You may not be able to do all of steps below  #
# in the same R session. I have had issues with #
# building the binary package after installing  #
# from source in the same session. YMMV.        #
#################################################



#Check the package (as for submission to CRAN).
devtools::check()


#############
# Make sure to increment the version number before running the following commands!

# Build source archive of package
devtools::build(packageLocation, path=paste0(baseDir, "packaged/") )

# Build binary archive of package
devtools::build(packageLocation, path=paste0(baseDir, "packaged/"), binary = TRUE)


#############
# Test installation from packaged files
install.packages( paste0(baseDir, "/packaged/CatContModel_0.9.0.tar.gz"), repos=NULL)
install.packages( paste0(baseDir, "/packaged/CatContModel_0.9.0.zip"), repos=NULL)


