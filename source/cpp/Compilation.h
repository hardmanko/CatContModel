#pragma once

#define COMPILING_WITH_RCPP 1

#ifdef COMPILING_WITH_RCPP

#include <Rcpp.h>

#define OVERRIDE

#else

#define COMPILING_WITH_CX 1

#include "CX.h"

#define OVERRIDE override

#endif