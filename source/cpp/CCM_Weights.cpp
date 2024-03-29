#include "CCM_Weights.h"

#include "CCM_ModelConfig.h"

//#include "CCM_Linear.h"
//#include "CCM_Circular.h"

// [[Rcpp::export]]
std::vector<double> CCM_CPP_CalculateWeights(std::string dataType, std::string weightType, double study, std::vector<double> activeCatMu, double catSelectivity) {

	CatCont::WeightsConfig wc;
	wc.dataType = CatCont::dataTypeFromString(dataType);
	wc.weightDist = CatCont::weightsDistributionFromString(weightType);
	wc.maxCategories = activeCatMu.size();

	CatCont::WeightsCalculator calc(wc);

	calc.calcWeights(study, activeCatMu, catSelectivity);

	return calc.copyFilledWeights();
}

namespace CatCont {

	WeightsDistribution weightsDistributionFromString(string weightsDistStr) {
		WeightsDistribution wd = WeightsDistribution::Default;
		if (weightsDistStr == "Nearest") {
			wd = WeightsDistribution::Nearest;
		}
		else if (weightsDistStr == "NearestInRange") {
			wd = WeightsDistribution::NearestInRange;
		}
#ifdef USING_CAT_SPLINE
		else if (weightsDistStr == "PlatSpline") {
			wd = WeightsDistribution::PlatSpline;
		}
#endif
		return wd;
	}

#ifdef USING_CAT_SPLINE
	LambdaVariant lambdaVariantFromString(string lambdaVariantStr) {

		LambdaVariant lambdaVariant = LambdaVariant::None;
		if (lambdaVariantStr == "CatWeightSum") {
			lambdaVariant = LambdaVariant::CatWeightSum;
		}
		else if (lambdaVariantStr == "InverseCatWeightSum") {
			lambdaVariant = LambdaVariant::InverseCatWeightSum;
		}
		else if (lambdaVariantStr == "Minus1to1") {
			lambdaVariant = LambdaVariant::Minus1to1;
		}
		return lambdaVariant;

	}
#endif



	///////////////////////
	// WeightsCalculator //
	///////////////////////

	WeightsCalculator::WeightsCalculator(const ModelConfiguration& modCfg) {
		WeightsConfig cfg;
		cfg.dataType = modCfg.dataType;
		cfg.maxCategories = modCfg.maxCategories;

		*this = WeightsCalculator(cfg);
	}

	WeightsCalculator::WeightsCalculator(const WeightsConfig& wConfig) :
		config(wConfig)
	{
		this->weights.resize(wConfig.maxCategories); // this always gives enough room
	}

	void WeightsCalculator::calcWeights(double study, const vector<double>& catMu, double catSelectivity) {

		this->weightCount = catMu.size();
		//this->weights.resize(this->weightCount); // See constructor

		if (catMu.size() == 0) {
			return;
		}

		switch (config.weightDist) {
		case WeightsDistribution::Default:
			_calcWeights_default_HVR17(study, catMu, catSelectivity);
			break;
		case WeightsDistribution::Nearest:
			_calcWeights_nearest(study, catMu);
			break;
		case WeightsDistribution::NearestInRange:
		{
			double catRange = catSelectivity;
			if (config.dataType == DataType::Circular) {
				double catRangeSD = CatCont::Circular::precRad_to_sdDeg(catSelectivity);
				catRange = CatCont::Circular::degreesToRadians(catRangeSD);
			}
			_calcWeights_nearestInRange(study, catMu, catRange);
		}
			break;
#ifdef USING_PLAT_SPLINE
		case WeightsDistribution::PlatSpline:
			_calcWeights_platSpline(study, catPar);
			break;
#endif
		}

	}

	double WeightsCalculator::sumWeights(void) const {
		double wSum = 0;
		for (size_t i = 0; i < this->weightCount; i++) {
			wSum += this->weights[i];
		}
		return wSum;
	}

	vector<double> WeightsCalculator::copyFilledWeights(void) const {
		return vector<double>(weights.begin(), weights.begin() + weightCount);
	}

	void WeightsCalculator::clearWeights(void) {
		this->weightCount = 0;
	}

	void WeightsCalculator::_calcWeights_default_HVR17(double study, const vector<double>& catMu, double catSelectivity) {

		if (catMu.size() == 1) {
			this->weights[0] = 1;
			return;
		}

		double wSum = 0;
		for (size_t i = 0; i < this->weightCount; i++) {
			double dens = 0;
			if (config.dataType == DataType::Linear) {
				//This distribution is not truncated: Category assignment is independent of the study/response space.
				dens = Linear::normalPDF(study, catMu[i], catSelectivity);
			}
			else if (config.dataType == DataType::Circular) {
				dens = vmLut.dVonMises(study, catMu[i], catSelectivity);
			}

			this->weights[i] = dens;

			wSum += dens;
		}

		_rescaleWeights(wSum);
		_zapSmallWeights();
	}

	void WeightsCalculator::_rescaleWeights(double wSum) {
		if (wSum < config.zeroSumCutoff) {

			// If wSum is tiny, give equal weights. This is rare.
			for (size_t i = 0; i < this->weightCount; i++) {
				this->weights[i] = 1.0 / this->weightCount;
			}
		} else {
			// Rescale weights to sum to 1
			for (size_t i = 0; i < this->weightCount; i++) {
				this->weights[i] /= wSum;
			}
		}
	}

	void WeightsCalculator::_zapSmallWeights(void) {
		if (config.zapsmallCutoff == 0) {
			return;
		}

		double wSum = 0;
		for (size_t i = 0; i < this->weightCount; i++) {
			if (this->weights[i] <= config.zapsmallCutoff) {
				this->weights[i] = 0;
			} else {
				wSum += this->weights[i];
			}
		}

		_rescaleWeights(wSum);
	}



	std::vector<double> WeightsCalculator::_studyMuDistance(double study, const vector<double>& catMu) const {
		std::vector<double> distances(catMu.size());

		for (size_t i = 0; i < catMu.size(); i++) {
			if (config.dataType == DataType::Linear) {
				distances[i] = abs(study - catMu[i]);
			}
			else if (config.dataType == DataType::Circular) {
				distances[i] = Circular::circularAbsoluteDistance(study, catMu[i]);
			}
		}

		return distances;
	}

	void WeightsCalculator::_calcWeights_nearest(double study, const vector<double>& catMu) {

		std::vector<double> distances = _studyMuDistance(study, catMu);

		size_t minDistI = 0;
		double minDist = std::numeric_limits<double>::max();
		for (size_t i = 0; i < catMu.size(); i++) {

			this->weights[i] = 0; // clear all weights

			if (distances[i] < minDist) {
				minDistI = i;
				minDist = distances[i];
			}
		}

		this->weights[minDistI] = 1;
	}

	void WeightsCalculator::_calcWeights_nearestInRange(double study, const vector<double>& catMu, double catRange) {
		std::vector<double> distances = _studyMuDistance(study, catMu);

		size_t minDistI = std::numeric_limits<size_t>::max();
		double minDist = std::numeric_limits<double>::max();
		for (size_t i = 0; i < catMu.size(); i++) {

			this->weights[i] = 0; // clear all weights

			if (distances[i] < catRange && distances[i] < minDist) {
				minDistI = i;
				minDist = distances[i];
			}
		}

		if (minDistI != std::numeric_limits<size_t>::max()) {
			this->weights[minDistI] = 1;
		}
		
	}

	/*
	void WeightsCalculator::_calcWeights_linear(double study, const vector<double>& catMu, double catSelectivity) {

		double wSum = 0;
		for (size_t i = 0; i < this->weightCount; i++) {
			//This distribution is not truncated: Category assignment is independent of the study/response space.
			double dens = Linear::normalPDF(study, catMu[i], catSelectivity);

			this->weights[i] = dens; // Set density for each.

			wSum += dens;
		}

		_rescaleWeights(wSum);
		_zapSmallWeights();
	}

	void WeightsCalculator::_calcWeights_circular(double study, const vector<double>& catMu, double catSelectivity) {

		double wSum = 0;
		for (size_t i = 0; i < this->weightCount; i++) {
			double dens = vmLut.dVonMises(study, catMu[i], catSelectivity);
			wSum += dens;
			this->weights[i] = dens;
		}

		_rescaleWeights(wSum);
		_zapSmallWeights();
	}
	*/

	namespace Circular {
		//OUT_weights must be an array of the correct size. 
		//This code is nasty for the sake of efficiency, avoiding lots of little allocations of weights vectors.
		//This function is called  in the likelihood function for every observation, so it gets a lot of use.
		void categoryWeights(double study, const MemFirst::CombinedParameters& par, double* OUT_weights) {

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
			}
			else {
				for (unsigned int i = 0; i < n; i++) {
					OUT_weights[i] /= densSum;
				}
			}
		}
	}

	namespace Linear {
		//OUT_weights must be an array of the correct size (par.cat.mu.size()). 
		//This code is nasty for the sake of efficiency, avoiding lots of little allocations of weights vectors.
		//This function is called  in the likelihood function for every observation, so it gets a lot of use.
		void categoryWeights(double study, const MemFirst::CombinedParameters& par, const LinearConfiguration& lc, double* OUT_weights) {

			unsigned int n = par.cat.mu.size();

			double densSum = 0;
			for (unsigned int i = 0; i < n; i++) {
				//This distribution is not truncated: Category assignment is independent of the study/response space.
				double d = normalPDF(study, par.cat.mu[i], par.cat.selectivity);
				densSum += d;
				OUT_weights[i] = d;
			}

			if (densSum < 1e-250) {

				//If densSum is tiny, give equal weights
				for (unsigned int i = 0; i < n; i++) {
					OUT_weights[i] = 1.0 / n;
				}
			}
			else {
				for (unsigned int i = 0; i < n; i++) {
					OUT_weights[i] /= densSum;
				}
			}
		}
	}


#ifdef USING_PLAT_SPLINE

	double WeightsCalculator::calcLambda(void) const {

		double wSum = this->sumWeights();

		double lambda = 0;
		switch (this->modCfg->lambdaVariant) {
		case LambdaVariant::None:
			break; // lambda = 0
		case LambdaVariant::CatWeightSum:
			lambda = clamp(wSum, 0, 1);
			break;
		case LambdaVariant::InverseCatWeightSum:
			lambda = 1 - clamp(wSum, 0, 1);
			break;
		case LambdaVariant::Minus1to1:
			lambda = clamp(wSum, 0, 1) * 2 - 1;
			break;
		}

		return lambda;
	}

	void WeightsCalculator::_calcWeights_platSpline(double study, const CategoryParameters& catPar) {
		// there may be linear and circular versions as well for this

		// get weights
		vector<double> ws = dPlatSplineWeights(study, catPar, *this->modCfg);

		/*
		Rcpp::Rcout << "_calcWeights_platSpline(): ws = " << std::endl;
		for (double w : ws) {
			Rcpp::Rcout << w << std::endl;
		}
		*/

		// copy weights (uggggg)
		//this->weightCount = ws.size();
		for (size_t i = 0; i < this->weightCount; i++) {
			this->weights[i] = ws[i];
		}

		/*
		Rcpp::Rcout << "_calcWeights_platSpline(): weights = " << std::endl;
		for (double w : this->weights) {
			Rcpp::Rcout << w << std::endl;
		}
		*/

		// No need to rescale weights for PlatSpline. In fact, 0 weight is ok.
		// Or maybe make sure total weight < 1, but I'm pretty sure that can't happen.
	}

#endif

#ifdef USING_PLAT_SPLINE

	/*
	double calcPlatSplineLambda(const std::vector<double>& weights, LambdaVariant variant) {

		double wSum = 0;
		for (const double& w : weights) {
			wSum += w;
		}

		double lambda = 0;
		switch (variant) {
		case LambdaVariant::None:
			break; // lambda = 0
		case LambdaVariant::CatWeightSum:
			lambda = clamp(wSum, 0, 1);
			break;
		case LambdaVariant::InverseCatWeightSum:
			lambda = 1 - clamp(wSum, 0, 1);
			break;
		case LambdaVariant::Minus1to1:
			lambda = clamp(wSum, 0, 1) * 2 - 1;
			break;
		}

		return lambda;
	}
	*/

	double dPlatSpline_distance(double x, double mu, bool abs, const ModelConfiguration& modCfg, bool degrees) {

		if (modCfg.dataType == DataType::Linear) {
			return std::abs(x - mu);
		}
		else {
			//d = Circular::circularAbsoluteDistance(x, mu, degrees);
			return Circular::circularDistance(x, mu, true, degrees);
		}

	}

	// why is this a lambda thing? why not just a function?
	std::function<double(double, double)> dPlatSpline_getDistanceFunction(const ModelConfiguration& modCfg) {
		// TODO: Optimize this, maybe by making distFun a class member

		std::function<double(double, double)> distFun; // Distance from x to mu

		if (modCfg.dataType == DataType::Linear) {
			distFun = [](double x, double mu) -> double {
				return x - mu;
			};
		}
		else if (modCfg.dataType == DataType::Circular) {
			distFun = [](double x, double mu) -> double {
				return Circular::circularSignedDistance(x, mu); // radians
			};
		}

		return distFun;
	}

	// Special cases handled:
	// 0 or 1 categories
	// 1 and 2+ zero-distance mus
	// All mu on the same side of x (linear data)
	//
	// TODO: Consider using a binning approach to finding indices. Basically convert x and mu to
	// ints (or whatever division into many steps). This makes more sense if study is a vector.
	std::vector<double> dPlatSplineWeights(double study, CategoryParameters catPar, const ModelConfiguration& modCfg) {

		if (modCfg.dataType == CatCont::DataType::Circular) {
			// clamp angles to [0, 2pi)
			study = Circular::clampAngle360(study);
			catPar.mu = Circular::clampAngle360(catPar.mu);
		}


		std::function<double(double, double)> distFun = dPlatSpline_getDistanceFunction(modCfg);

		auto singleMuDensity =
			[&distFun](double study, const CategoryParameters& catPar, size_t muIndex) -> double
		{

			double d = distFun(study, catPar.mu[muIndex]);
			//Rcpp::Rcout << "singleMuDensity(): d = " << d << std::endl;
			d = std::abs(d);

			double platHW = catPar.platHW[muIndex];
			double splineHW = catPar.selectivity;

			double rval = dPlatSpline(d, platHW, splineHW);
			//Rcpp::Rcout << "singleMuDensity(): rval = " << rval << std::endl;
			return rval;
		};

		auto approxZero = [](double dist) -> bool {
			return std::abs(dist) < 1e-12;
		};



		// rval is of correct dimension and initialized to 0s.
		std::vector<double> rval(catPar.nCat(), 0.0); // Default 0 density

		/*
		auto printRval = [&rval]() {
			Rcpp::Rcout << "weights: " << std::endl;
			for (double w : rval) {
				Rcpp::Rcout << w << std::endl;
			}
			Rcpp::Rcout << std::endl;
		};
		*/

		// Special cases for 0 or 1 categories.
		if (catPar.nCat() == 0) {
			return rval;
			//return std::vector<double>();
		}
		else if (catPar.nCat() == 1) {

			rval.front() = singleMuDensity(study, catPar, 0);
			//printRval();
			return rval;
		}

		////////
		// Calculate distances and look for 0 distance

		std::vector<double> singedDists(catPar.nCat());
		std::vector<size_t> zeroIndices;

		for (size_t j = 0; j < catPar.mu.size(); j++) {
			singedDists[j] = distFun(study, catPar.mu[j]);

			if (approxZero(singedDists[j])) {
				zeroIndices.push_back(j);
			}
		}

		/*
		Rcpp::Rcout << "SignedDists:" << std::endl;
		for (size_t i = 0; i < singedDists.size(); i++) {
			Rcpp::Rcout << singedDists[i] << std::endl;
		}

		Rcpp::Rcout << "zeroIndices:" << std::endl;
		for (size_t i = 0; i < zeroIndices.size(); i++) {
			Rcpp::Rcout << zeroIndices[i] << std::endl;
		}
		*/

		if (zeroIndices.size() >= 2) {
			//Rcpp::Rcout << "2+ zero indices: Setting the zero distance cats to equal weight." << std::endl;
			for (size_t ind : zeroIndices) {
				rval[ind] = 1.0 / zeroIndices.size();
			}
			return rval;
		}



		//////
		// Get indices of 2 nearest mus

		double lowerDist = std::numeric_limits<double>::lowest();
		double upperDist = std::numeric_limits<double>::max();
		size_t lowerIndex = std::numeric_limits<size_t>::max();
		size_t upperIndex = std::numeric_limits<size_t>::max();

		if (zeroIndices.size() == 1) {

			//Rcpp::Rcout << "1 zero distance mu" << std::endl;

			// One of the distances is 0. Find the other nearest mu.

			double nearestDist = std::numeric_limits<double>::max();
			size_t nearestIndex = std::numeric_limits<size_t>::max();

			for (size_t j = 0; j < catPar.mu.size(); j++) {

				if (j == zeroIndices.front()) {
					continue;
				}

				if (std::abs(singedDists[j]) < std::abs(nearestDist)) {
					nearestIndex = j;
					nearestDist = singedDists[j];
				}
			}

			// Special case: If one of the distances is 0,
			// then the other side is the closest other mu.
			if (nearestDist < 0) {
				// If the nearest non-zero distance is negative, so lower is nearest and upper is zero
				lowerDist = nearestDist;
				lowerIndex = nearestIndex;

				upperDist = 0.0;
				upperIndex = zeroIndices.front();
			}
			else {
				// Nearest positive, so lower is zero and upper is nearest
				lowerDist = 0.0;
				lowerIndex = zeroIndices.front();

				upperDist = nearestDist;
				upperIndex = nearestIndex;
			}

		}
		else {
			// None of the distances are approx 0. Find the nearest indices normally.
			//Rcpp::Rcout << "0 zero distance mus" << std::endl;

			for (size_t j = 0; j < catPar.nCat(); j++) {
				double dist = singedDists[j];

				// If distance is negative and less negative than nearest, save it
				if (dist < 0.0 && dist > lowerDist) {
					lowerDist = dist;
					lowerIndex = j;
				}
				else if (dist > 0.0 && dist < upperDist) {
					upperDist = dist;
					upperIndex = j;
				}
			}

		}

		/*
		Rcpp::Rcout << "lowerDist: " << lowerDist << std::endl;
		Rcpp::Rcout << "lowerIndex: " << lowerIndex << std::endl;
		Rcpp::Rcout << "upperDist: " << upperDist << std::endl;
		Rcpp::Rcout << "upperIndex: " << upperIndex << std::endl;
		*/

		// This is impossible if nCat() > 0, which it is, so this is an error case.
		//if (lowerIndex == std::numeric_limits<size_t>::max() && upperIndex == std::numeric_limits<size_t>::max()) {
		//	return rval;
		//}

		// For linear data, you can have all categories on only 1 side of x,
		// so check for those two cases. If all mus are on the same side of x,
		// only get a density for the nearest mu.
		if (lowerIndex == std::numeric_limits<size_t>::max()) {
			//Rcpp::Rcout << "No lower index found" << std::endl;
			rval[upperIndex] = singleMuDensity(study, catPar, upperIndex);
			return rval;
		}
		if (upperIndex == std::numeric_limits<size_t>::max()) {
			//Rcpp::Rcout << "No upper index found" << std::endl;
			rval[lowerIndex] = singleMuDensity(study, catPar, lowerIndex);
			return rval;
		}

		//assert(upperDist >= 0 && lowerDist <= 0);
		if (!(upperDist >= 0 && lowerDist <= 0)) {
			Rcpp::Rcout << "Upper and lower distances are not on their correct side." << std::endl;
		}

		// Copy out parameters of interest (so they can be modified).
		// PlatWidth, PW, and SplineWidth, SW, are both half widths.
		double lowerPW = catPar.platHW[lowerIndex];
		double upperPW = catPar.platHW[upperIndex];
		double lowerSW = catPar.selectivity; // [lowerIndex]?
		double upperSW = catPar.selectivity;

		double lowerPSW = lowerPW + lowerSW;
		double upperPSW = upperPW + upperSW;

		double totalPSW = lowerPSW + upperPSW;

		double maxPSW = upperDist - lowerDist;

		// Modifies its arguments
		auto copeWithLoss = [](double loss, double& platHW, double& splineHW) -> bool {

			//assert(loss > 0 && loss <= platHW + splineHW);

			// If the loss is less than the plateau, just shorten the plateau.
			if (loss <= platHW) {
				platHW -= loss;
				return true;
			}

			// If the loss is greater than platHW, reduce the loss by platHW
			loss -= platHW;
			// reduce platHW to 0
			platHW = 0;
			// and reduce splineHW by the remaining loss
			splineHW -= loss;
			// It is not possible for the remaining loss to be greater than splineHW if 0 <= loss <= (platHW + splineHW)

			return true;
		};

		if (totalPSW > maxPSW) {
			//Rcpp::Rcout << "Coping with loss. totalPSW = " << totalPSW << ", maxPSW = " << maxPSW << std::endl;

			// rescale parameters
			double newScale = maxPSW / totalPSW; // totalPSW > maxPSW  <=>  newScale < 1
			//double lossP = 1 - maxPSW / totalPSW;

			// loss is bounded between 0 and lowerPSW
			double lowerLoss = lowerPSW * (1 - newScale);
			double upperLoss = upperPSW * (1 - newScale);

			// Note reference parameters modified in function
			copeWithLoss(lowerLoss, lowerPW, lowerSW);
			copeWithLoss(upperLoss, upperPW, upperSW);
		}

		rval[lowerIndex] = dPlatSpline(std::abs(lowerDist), lowerPW, lowerSW);

		rval[upperIndex] = dPlatSpline(upperDist, upperPW, upperSW);
		//rval[upperIndex] = 1 - rval[lowerIndex]; // This works only if the distributions are pushed up against one another

		return rval;
	}


	// z <= 0 is the high end of the spline, z == 1 is the middle, and z >= 2 is the low end
	double zeroDerivativeCubicSplineDensity(double z) {


		/*
		if (z < 0.0) {
		//Enforce z >= 0?
			return NA_REAL;
		} else if (z > 2.0) {
			return 0.0; // Off the end of the tail
		}
		*/

		if (z < 0.0) {
			// On the plateau, probably
			return 1.0;
		}
		else if (z > 2.0) {
			// Off the end of the tail
			return 0.0;
		}

		// These coefficients create a spline with slope 0 at either end.
		// I used the method described here: https://www.tf.uni-kiel.de/matwis/amat/comp_math/kap_1/backbone/r_se16.html
		// See private/AltWeightsFunction_Spline.R for the specific code I used.
		//static const double coef[4] = { 1.0, 0.0, -0.75, 0.25 };
		//
		//double dens = 0;
		//for (size_t i = 0; i < 4; i++) {
		//	dens += pow(z, (double)i) * coef[i];
		//}

		// The above loop unrolled (half as many pow calls)
		double dens = 1.0; // 1 * z^0 = 1
		//dens += coef[1] * z; // 0 * z^1 = 0
		dens += -0.75 * pow(z, 2.0);
		dens += 0.25 * pow(z, 3.0);

		return dens;
	}

	// distFromPlatEdge has already had the plateau subtracted off
	double zeroDerivativeCubicSplineDensity2(double distFromPlatEdge, double splineHW) {
		double z = distFromPlatEdge / splineHW;
		return zeroDerivativeCubicSplineDensity(z);
	}

	// absDist is distance from center of plateau (i.e. mu) to the study value.
	// platHW is the plateau half width
	// splineHW is the spline half width ("standard deviation")
	double dPlatSpline(double absDist, double platHW, double splineHW) {

		if (absDist < platHW) {
			return 1.0; // On the plateau
		}
		else if (absDist > platHW + 2 * splineHW) {
			return 0.0; // Off the tail of the spline
		}

		double remD = absDist - platHW;

		// TODO: Decide where to check this, in this function or in zeroDerivativeCubicSplineDensity
		/*
		if (remD < 0) {
			return 1.0; // On the plateau
		} else if (remD >= 2 * splineHW) {
			return 0.0; // Off the tail of the spline
		}
		*/

		return zeroDerivativeCubicSplineDensity2(remD, splineHW);
	}

	double dPlatSplineFull(double x, double mu, double platHW, double splineHW, bool linear, bool degrees) {

		double d;
		if (linear) {
			d = std::abs(x - mu);
		}
		else {
			//d = Circular::circularAbsoluteDistance(x, mu, degrees);
			d = Circular::circularDistance(x, mu, true, degrees);
		}

		return dPlatSpline(d, platHW, splineHW);
	}

#endif

} // namespace CatCont