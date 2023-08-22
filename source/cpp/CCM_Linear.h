#pragma once

#include "CCM_Main.h"

namespace CatCont {

	namespace Linear {

		double normalPDF(double x, double mu, double sd);
		double normalCDF(double x, double mu, double sd);
		double dtnorm_denominator(double mu, double sd, double lower, double upper);
		double dtnorm(double x, double mu, double sd, double lower, double upper, bool log_ = false);
		vector<double> dtnorm(const vector<double>& x, double mu, double sd, double lower, double upper, bool log_);
		double dtnorm_noBoundCheck(double x, double mu, double sd, double lower, double upper);

		// For within-item
		double combineSDs(double pContWithin, double contSD, double catSD);


	}
}