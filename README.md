# Introduction

This R package implements the models used by Hardman, Vergauwe, and Ricker (2016). These models are difficult to implement and difficult to use. This package simplifies the process of using the models.

This package is under development, as are the models. If you run into bugs or unusual model behavior, please report them in the issue tracker (click on the "Issues" tab).


# Installation

I would recommend that you download a release of the package (see the releases tab). A release contains the whole repository at a particular point in time, which includes the introduction documentation pdf and examples. The release itself is not an installable R package, but it contains, in the "packaged" directory, compiled versions of the package that can be installed by R.

Once you have downloaded a release, unzip it somewhere on your computer. From there, the easiest way to install the package is to use the "installPackage.R" script in the root directory. Once you set a path in it to point to one of the installable packages, it installs a version of the package that should work for most people. If you are on Linux or OSx, you should see the Installing from Source section below.

Alternately, if you don't want the introduction documentation or examples (you will still get function documentation), you can just download a compiled version of the package from the "packaged" subdirectory of this repository. A .zip extension on a package file there indicates that it is a binary (compiled) version of the package, which is what Windows users will want, whereas a .tar.gz indicates a source version of the package, which is what Linux and OSx users will want.

If you need to install from source, such as if you are on Linux or OSx, see the "Installing from Source" section below.


## Installing from Source

Source versions of the package can be found in the "packaged" subdirectory and have a ".tar.gz" file extension. If you are on Linux of OSx, you may need to install the package from source as don't have compiled versions of the package for those platforms.

Basically, installing from source requires that you have a C++ compiler that R knows how to work with configured on your computer. Below are some instructions to help with it.

You will need to have a C++ compiler WITH C++11 SUPPORT installed on your computer. This is not something that comes with R.
- On Windows, you can download RTools (not an R package), which will install a C++ compiler. See https://cran.r-project.org/bin/windows/Rtools/
- On OSx, read this: http://seananderson.ca/2013/11/18/rcpp-mavericks.html
- On linux, you probably already have a C++ compiler installed.

If the compiler you have is old (or cranky), you may need to explicitly tell R to compile C++ code with C++11 support. 
To enable C++11 support, you will need a Makevars file with one of the following lines in it.

If you have an older compiler:

CXXFLAGS+=-std=c++0x

If you have a newer compiler:

CXXFLAGS+=-std=c++11

If you are uncertain, you may need to try both. The compiler that comes with older versions of RTools is really old, so it needs "c++0x". RTools for R 3.3 can use "c++11".

The location of the Makevars file and even what name it should have is a little bit of a mystery.
See https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Using-Makevars
That page supposedly explains what Makevars is, where it should go, and what name it has,
but I disagree with the "experts" that the page actually does any of those things.
I had to make a few educated guesses to get it to work.

At this point, you should be able to install from source. Make sure to install the dependencies using the "InstallPackage.R" script before installing this package.


# Usage

Start by reading the introduction to use of this package in docs/introduction/Introduction.pdf.

Some examples that use simulated data can be found in "examples". They are not necessarily supposed to be run in a totally linear fashion and they are used for testing the package, so they are not the cleanest, but they do give some usage examples.


# Building the Package

This repo contains packaged versions of the package that can be installed by R fairly directly. If you just want to use the package as is, you just need to install those packaged versions (see Installation).

If you want to develop/modify the package, you will need to rebuild the package after making any changes. This repository is not itself in the correct format for an R package, but rather a collection of files that are combined together into an R package with buildPackage.R, which is in the root directory. You can typically just run most of that file, but I have had some problems with the last couple of steps all working in the same R session, so read the comments toward the end.


# License

This package is released under the MIT license. See LICENSE.md for more information.


