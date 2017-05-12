# Introduction

This R package implements the models used by Hardman, Vergauwe, and Ricker (2016) [manuscript link](http://kylehardman.com/Content/StaticPages/Publications/Files/Hardman%20Vergauwe%20and%20Ricker%20(2017)%20-%20Manuscript.pdf). These models are difficult to implement and difficult to use. This package simplifies the process of using the models.

This package is under development, as are the models. If you run into bugs or unusual model behavior, please report them in the issue tracker (click on the "Issues" tab).


# Installation

To install the package, you will need two things:

1. An installable copy of the package.
2. The file "installPackage.R", which can be found in the same directory as this readme file.

These are discussed in the next two sections.

## Installable copy of the package

This repository is not itself installable by R, as many R package repositories are. As a result, you can't use, e.g. `devtools::install_github` to install it.

You can get an installable copy of the package from the releases of the package, see the [releases tab](https://github.com/hardmanko/CatContModel/releases). For versions of the package older than 0.7.4, installable copies of the package can be found in the "packaged" directory of this repository.

Installable package files have filenames beginning with `CatContModel_V.V.V`, where `V.V.V` is the version number. 

The type of installable package that you install depends on your operating system. For Windows and OSx, binary versions of the package that do not require C++ compilation are available and highly recommended.

Different package types are identified by their file extension:

+ Windows binary: `.zip` file extension.
+ OSx binary: `.tgz` file extension.
+ Source: `tar.gz` file extension. Can be installed on any operating system, but requires a properly-configured C++ compiler. See the Installing from Source section below.


## installPackage.R

Once you have an installable version of the package, the file [installPackage.R](https://github.com/hardmanko/CatContModel/blob/master/installPackage.R) contains code for installing this package. It installs package dependencies and gives information on installing this package. Comments in that file indicate how to use it. 

# Usage

## Manual

Once you have finished installing the package, read the introduction to use of this package. The latest version of the manual can be found in [docs/introduction/Introduction.pdf](https://github.com/hardmanko/CatContModel/raw/master/docs/introduction/Introduction.pdf).

The manual should also be installed as a vignette when you install the package. Open it from R with.

```{r}
vignette("Introduction", "CatContModel")
```

## Examples

Some examples that use simulated data can be found in the [examples](https://github.com/hardmanko/CatContModel/tree/master/examples) directory. They are not necessarily supposed to be run in a totally linear fashion and they are used for testing the package, so they are not the cleanest, but they do give some usage examples with comments.

The easiest way to get the examples is to download the whole repository and unzip it somewhere on your computer. See the green Clone or Download button in the top right.


# Installing from Source

Source versions of the package can be found in the "packaged" subdirectory and have a ".tar.gz" file extension. If you are on Linux or an unusual OS (like Solaris), you need to install the package from source as I don't have compiled versions of the package for those platforms.

Installing from source requires that you have a C++ compiler that R knows how to work with configured on your computer. Below are some instructions to help with it.

You will need to have a C++ compiler WITH C++11 SUPPORT installed on your computer. This is not something that comes with R.
- On Windows, you can download RTools (not an R package), which will install a C++ compiler. See https://cran.r-project.org/bin/windows/Rtools/
- On OSx, read this: http://seananderson.ca/2013/11/18/rcpp-mavericks.html
- On Linux, you may already have an appropriate C++ compiler installed.

If the compiler you have is old (or cranky), you may need to explicitly tell R to compile C++ code with C++11 support. 
To enable C++11 support, you will need a Makevars file with one of the following lines in it.

If you have an older compiler:

    CXXFLAGS+=-std=c++0x

If you have a newer compiler:

    CXXFLAGS+=-std=c++11

If you are uncertain, you may need to try both. The compiler that comes with older versions of RTools is really old, so it needs "c++0x". RTools (Windows) for R 3.3 can use "c++11".

The location of the Makevars file and even what name it should have is a little bit of a mystery.
See https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Using-Makevars
That page supposedly explains what Makevars is, where it should go, and what name it has,
but I disagree with the "experts" that the page actually does any of those things.
I had to make a few educated guesses to get it to work.

At this point, you should be able to install from source. Make sure to install the dependencies using the "InstallPackage.R" script before installing this package.


# Building the Package

This repo contains packaged versions of the package that can be installed by R directly. If you just want to use the package as is, you just need to install those packaged versions.

If you want to develop/modify the package, you will need to rebuild the package after making any changes. This repository is not itself in the correct format for an R package, but rather a collection of files that are combined together into an R package with buildPackage.R, which is in the root directory. You can typically just run most of that file, but I have had some problems with the last couple of steps all working in the same R session, so read the comments toward the end.


# License

This package is released under the MIT license. See LICENSE.md for more information.


