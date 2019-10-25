#pragma once

#include "Compilation.h"

#ifdef COMPILING_WITH_CX
#include "ofxRmath.h"
#endif

#include <vector>
#include <cmath>

#if COMPILING_WITH_RCPP
#define GS_COUT Rcpp::Rcout
#endif

#if COMPILING_WITH_CX
#define GS_COUT std::cout
#endif


double normal_varPostSample(const std::vector<double>& y, double mu, double a0, double b0);

double normal_muPostSample(double mean, double N, double data_var, double mu0, double var0);
double normal_muPostSample(const std::vector<double>& y, double data_var, double mu0, double var0);

//This has not been properly tested!
double inverseGammaLogLikelihood(double x, double shape, double scale);

double binomialProportionalLogLikelihood(double successes, double trials, double p);

double multinomialProportionalLogLikelihood(const std::vector<unsigned int>& x, const std::vector<double>& p);

double logitTransform(double unitInterval);
double logitInverse(double fullScale);

#ifdef COMPILING_WITH_CX
void outputDataVector(std::string filename, std::string header, std::vector<double> values);
#endif