# Introduction

This R package implements the models used by Hardman, Vergauwe, and Ricker (2016). These models are difficult to implement and difficult to use. This package simplifies the process of using the models.

This package is under development, as are the models. If you run into bugs or unusual model behavior, please report them in the issue tracker (click on the "Issues" tab).


# Installation

I would recommend that you download a release of the package (see the releases tab). A release contains the whole repository at a particular point in time, which includes the introduction documentation pdf and examples. The release itself is not an installable R package, but it contains, in the "packaged" directory, compiled versions of the package that can be installed by R.

Once you have downloaded a release, unzip it somewhere on your computer. From there, the easiest way to install the package is to use the "installPackage.R" script in the root directory. Once you set a path in it to point to one of the installable packages, it installs a version of the package that should work for most people. 

Alternately, if you don't want the introduction documentation or examples (you will still get function documentation), you can just download a compiled version of the package from the "packaged" subdirectory of this repository. A .zip extension on a packaged file there indicates that it is a binary (compiled) version of the package, which is almost certaintly what you want, not the .tar.gz, which is a source version of the package.  You should still use "InstallPackage.R" if you go this route.

If you need to install from source, see the "Installing from Source" section below.


# Usage

Start by reading the introduction to use of this package in docs/introduction/Introduction.pdf.

Some examples that use simulated data can be found in "examples". They are not necessarily supposed to be run in a totally linear fashion and they are used for testing the package, so they are not the cleanest, but they do give some usage examples.

Simulated sample data and analyses on those data sets can be found in the example folders.


# Building the Package

This repo contains packaged versions of the package that can be installed by R fairly directly. If you just want to use the package as is, you just need to install those packaged versions (see Installation).

If you want to develop/modify the package, you will need to rebuild the package after making any changes. This repository is not itself in the correct format for an R package, but rather a collection of files that are combined together into an R package with buildPackage.R, which is in the root directory. You can typically just run that whole file, but sometimes the


# Installing from Source

Source versions of the package can be found in the "packaged" subdirectory.

Basically, installing from source requires that you have a C++ compiler that R knows how to work with configured on your computer. Below are some instructions to help with it.

You will need to have a C++ compiler WITH C++11 SUPPORT installed on your computer. This is not something that comes with R.
- On Windows, you can download RTools (not a package), which will install a C++ compiler. See https://cran.r-project.org/bin/windows/Rtools/
- On OSx, read this: http://seananderson.ca/2013/11/18/rcpp-mavericks.html
- On linux, you probably already have a C++ compiler installed.

If the compiler you have is old, like the compiler that comes with RTools, 
you may need to explicitly tell R to compile C++ code with C++11 support. 
If using an old version of GCC, you will need a Makevars file with the following line in it:

CXXFLAGS+=-std=c++0x

OR if your compiler is new enough, it will be

CXXFLAGS+=-std=c++11

The compiler that comes with older versions of RTools is really old, so it needs "c++0x". RTools for
R 3.3 can use "c++11".

The location of the Makevars file and even what name it should have is a little bit of a mystery.
See https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Using-Makevars
That page supposedly explains what Makevars is, where it should go, and what name it has,
but I disagree with the "experts" that the page actually does that. I had to make a few educated
guesses to get it to work.

At this point, you should be able to install from source. Make sure to install the dependencies using the "InstallPackage.R" script.


# License

This package is released under the MIT license. See LICENSE.md for more information.


# Versions

0.6.1 - Minor bugfix in a plotting function. Additional documentation for the linear dataType.

0.6.0 - Added linear data type. Additional changes.

0.5.0 - Initial release.
