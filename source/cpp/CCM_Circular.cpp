#include "CCM_Circular.h"

namespace CatCont {
namespace Circular {

double degreesToRadians(double degrees) {
	return degrees * PI / 180.0;
}
double radiansToDegrees(double radians) {
	return radians * 180.0 / PI;
}

vector<double> degreesToRadians(const vector<double>& degrees) {
	vector<double> radians(degrees.size());
	for (unsigned int i = 0; i < degrees.size(); i++) {
		radians[i] = degrees[i] * PI / 180.0;
	}
	return radians;
}

vector<double> radiansToDegrees(const vector<double>& radians) {
	vector<double> degrees(radians.size());
	for (unsigned int i = 0; i < radians.size(); i++) {
		degrees[i] = radians[i] * 180.0 / PI;
	}
	return degrees;
}

double precRad_to_sdDeg(double precRad) {
	double sdRad = 1.0 / sqrt(precRad); //precision to sd
	return sdRad * 180.0 / PI; //rad to deg
}

double sdDeg_to_precRad(double sdDeg) {
	double sdRad = sdDeg * PI / 180.0; //degrees to rad
	return 1.0 / (sdRad * sdRad); //sd to precision
}



double circularMean(double rad1, double rad2) {
	vector<double> rads = { rad1, rad2 };
	//rads[0] = rad1;
	//rads[1] = rad2;
	return circularMean(rads);
}

double circularMean(const vector<double>& radians) {
	vector<double> weights(radians.size(), 1.0 / radians.size());
	return circularMean(radians, weights);
}

// weights should sum to 1
double circularMean(const vector<double>& radians, const vector<double>& weights) {
	double cosx = 0;
	double sinx = 0;
	for (size_t i = 0; i < radians.size(); i++) {
		cosx += cos(radians[i]) * weights[i];
		sinx += sin(radians[i]) * weights[i];
	}
	double atn = atan2(sinx, cosx);
	if (atn < 0) {
		atn += 2 * PI;
	}
	return atn;
}


// 0 <= rval < 2pi
double clampAngle360(double x) {

	x = std::fmod(x, 2*PI);

	if (x < 0) {
		x += 2*PI;
	}

	return x;
}

// -pi <= rval < pi
double clampAngle180(double x) {

	x = std::fmod(x, 2 * PI);

	if (x < -PI) {
		x += 2 * PI;
	} else if (x >= PI) {
		x -= 2 * PI;
	}

	return x;
}

vector<double> clampAngle180(const vector<double>& xs) {
	vector<double> rval(xs.size());
	for (size_t i = 0; i < xs.size(); i++) {
		rval[i] = clampAngle180(xs[i]);
	}
	return rval;
}

vector<double> clampAngle360(const vector<double>& xs) {
	vector<double> rval(xs.size());
	for (size_t i = 0; i < xs.size(); i++) {
		rval[i] = clampAngle360(xs[i]);
	}
	return rval;
}

double clampAngle(double x, bool pm180, bool degrees) {

	if (degrees) {
		x = degreesToRadians(x);
	}

	double rval = pm180 ? clampAngle180(x) : clampAngle360(x);

	if (degrees) {
		rval = radiansToDegrees(rval);
	}

	return rval;

	/*

	double fcirc = degrees ? 360 : 2 * PI;
	double hcirc = degrees ? 180 : PI;

	x = std::fmod(x, fcirc);

	if (pm180) {
		// -180 <= x <= 180
		if (x <= -hcirc) {
			x += fcirc;
		} else if (x > hcirc) {
			x -= fcirc;
		}
	} else {
		// 0 <= x <= 360
		if (x < 0) {
			x += fcirc;
		}
	}

	return x;
	*/
}


// assume radians
double circularSignedDistance(double x, double y) {
	return clampAngle180(y - x);
}

double circularAbsoluteDistance(double x, double y) {
	return std::abs(circularSignedDistance(x, y));
}

double circularDistance(double x, double y, bool absDist, bool degrees) {
	if (degrees) {
		x = degreesToRadians(x);
		y = degreesToRadians(y);
	}

	double d = circularSignedDistance(x, y);
	if (absDist) {
		d = std::abs(d);
	}

	if (degrees) {
		d = radiansToDegrees(d);
	}

	return d;
}



//OUT_weights must be an array of the correct size. 
//This code is nasty for the sake of efficiency, avoiding lots of little allocations of weights vectors.
//This function is called  in the likelihood function for every observation, so it gets a lot of use.
void categoryWeights(double study, const CombinedParameters& par, double* OUT_weights) {

	unsigned int n = par.cat.mu.size();

	double densSum = 0;
	for (unsigned int i = 0; i < n; i++) {
		double d = vmLut.dVonMises(study, par.cat.mu[i], par.cat.selectivity);
		densSum += d;
		OUT_weights[i] = d;
	}

	if (densSum < 1e-250) {

		//If densSum is tiny, give equal weights
		for (unsigned int i = 0; i < n; i++) {
			OUT_weights[i] = 1.0 / n;
		}
	} else {
		for (unsigned int i = 0; i < n; i++) {
			OUT_weights[i] /= densSum;
		}
	}
}

// See Equation 25 in the Appendix. This is a rewritten form. sigma^2 = 1 / kappa
double combineKappas(double contKappa, double catKappa, double pCont) {

	// This is based on the idea that the study and category locations both have some amount of
	// noise added to them, then are averaged. Thus, the average has reduced variance.
	double combinedVariance = (pow(pCont, 2) / contKappa) + (pow((1 - pCont), 2) / catKappa);

	double newKappa = 1 / combinedVariance; //convert variance back to precision

	// Make sure the newKappa is in range. This is invisible to the rest of the system!
	// However, it is also very rare that this is needed, at least for my data.
	newKappa = std::min(newKappa, vmLut.getMaxKappa());

	return newKappa;
}


vector<double> zlLikelihood(const zlParameters& par, const ConditionData& data) {

	const double unifDens = 1 / (2 * PI);
	const double guessDens = (1 - par.pMem) * unifDens;

	vector<double> likelihoods(data.study.size());

	for (unsigned int i = 0; i < data.study.size(); i++) {
		double memDens = par.pMem * vmLut.dVonMises(data.response[i], data.study[i], par.contSD);

		likelihoods[i] = memDens + guessDens;

		likelihoods[i] *= (PI / 180.0);
	}
	return likelihoods;
}

vector<double> betweenAndWithinLikelihood(const CombinedParameters& par, const ConditionData& data, ModelVariant modelVariant) {

	if (modelVariant == ModelVariant::ZL) {

		zlParameters zlPar;
		zlPar.pMem = par.pMem;
		zlPar.contSD = par.contSD;

		return zlLikelihood(zlPar, data);
	}

	bool calculateWithinComponent = (modelVariant == ModelVariant::BetweenAndWithin) || (modelVariant == ModelVariant::WithinItem);
	bool calculateBetweenComponent = (modelVariant == ModelVariant::BetweenAndWithin) || (modelVariant == ModelVariant::BetweenItem);

	double pBetween = par.pBetween;
	if (modelVariant == ModelVariant::BetweenItem) {
		pBetween = 1;
	} else if (modelVariant == ModelVariant::WithinItem) {
		pBetween = 0;
	}

	const double unifDens = 1 / (2 * PI);

	vector<double> weights(par.cat.mu.size());

	//Precalculate this because it doesn't depend on the data
	// This variable naming SUCKS. These are not SDS, they have been converted to kappas.
	double combinedWithinKappa = Circular::combineKappas(par.contSD, par.cat.SD, par.pContWithin);


	//When calculating the within component predicted response center, weight the 
	//categorical and continuous components of the response by pContWithin
	vector<double> cmLoc(2);
	vector<double> cmWeights(2);
	cmWeights[0] = par.pContWithin;
	cmWeights[1] = 1 - par.pContWithin;


	vector<double> likelihoods(data.study.size());

	for (unsigned int i = 0; i < data.study.size(); i++) {

		//Category weights apply to both the between and with components of the model
		categoryWeights(data.study[i], par, weights.data());

		//Within component
		double withinDensity = 0;

		if (calculateWithinComponent) {

			if (par.cat.mu.size() == 0) {

				withinDensity = vmLut.dVonMises(data.response[i], data.study[i], par.contSD);

			} else {

				cmLoc[0] = data.study[i]; //The continuous component of the predicted response center

				for (unsigned int j = 0; j < weights.size(); j++) {

					cmLoc[1] = par.cat.mu[j]; //The categorical component of the center

					double predictedCenter = circularMean(cmLoc, cmWeights);

					double thisDens = vmLut.dVonMises(data.response[i], predictedCenter, combinedWithinKappa);

					withinDensity += weights[j] * thisDens;
				}

			}
		}


		double betweenCatDensity = 0;
		double catGuessDens = 0;

		//For both between and guessing. This is always calculated for guessing even if between isn't calculated.
		for (unsigned int j = 0; j < weights.size(); j++) {
			double dens = vmLut.dVonMises(data.response[i], par.cat.mu[j], par.cat.SD);

			betweenCatDensity += weights[j] * dens;
			catGuessDens += dens;
		}

		//Between component
		double betweenDensity = 0;

		if (calculateBetweenComponent) {
			double contDensity = vmLut.dVonMises(data.response[i], data.study[i], par.contSD);
			betweenDensity = (par.pContBetween * contDensity) + ((1 - par.pContBetween) * betweenCatDensity);
		}


		//Guessing component (sort of). catGuessDens is calculated above.
		if (par.cat.mu.size() > 0) {
			catGuessDens /= par.cat.mu.size();
		}

		double guessingDensity = par.pCatGuess * catGuessDens + (1 - par.pCatGuess) * unifDens;


		//Combine components
		double memoryDensity = pBetween * betweenDensity + (1 - pBetween) * withinDensity;

		double likelihood = par.pMem * memoryDensity + (1 - par.pMem) * guessingDensity;

		//The PI / 180 scales the likelihood so that it would be correct if degrees had been used instead of radians
		likelihoods[i] = likelihood * (PI / 180.0);

	}

	return likelihoods;

}

double betweenAndWithinLL(const CombinedParameters& par, const ConditionData& data, ModelVariant modelVariant) {
	vector<double> likelihoods = betweenAndWithinLikelihood(par, data, modelVariant);

	double llSum = 0;
	for (unsigned int i = 0; i < likelihoods.size(); i++) {
		llSum += log(likelihoods[i]);
	}

	return llSum;
}

double betweenAndWithinNll(const CombinedParameters& par, const ConditionData& data, ModelVariant modelVariant) {

	return -betweenAndWithinLL(par, data, modelVariant);

}

} // namespace Circular
} // namespace CatCont