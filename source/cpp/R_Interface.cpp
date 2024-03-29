#include "R_Interface.h"

#include "MF_RunModel.h"

#ifdef USING_CAT_SPLINE
#include "CCM_Weights.h"
#endif

#ifdef COMPILING_WITH_RCPP

vector<ParameterList> transposePosteriors(Rcpp::List posteriors) {

	vector<ParameterList> postIter;

	vector<string> postNames = posteriors.names();

	if (postNames.size() == 0) {
		return postIter;
	}

	size_t totalIterations = convertNumericVector(posteriors[postNames[0]]).size();
	postIter.resize(totalIterations);

	for (size_t i = 0; i < postNames.size(); i++) {

		string n = postNames[i];
		Rcpp::NumericVector thisPost = posteriors[n];

		for (size_t j = 0; j < totalIterations; j++) {
			postIter[j][n] = thisPost[j];
		}
	}

	return postIter;
}

std::vector<std::string> convertCharVector(const Rcpp::CharacterVector& cv) {
	return std::vector<std::string>(cv.begin(), cv.end());
}

std::vector<double> convertNumericVector(const Rcpp::NumericVector& nv) {
	return std::vector<double>(nv.begin(), nv.end());
}



// [[Rcpp::export]]
Rcpp::List CCM_CPP_runParameterEstimation(
	Rcpp::DataFrame data,
	Rcpp::List modelConfig,
	Rcpp::List runConfig,
	Rcpp::List equalityConstraints
)
{
	// This is for MemFirst only
	CatCont::MemFirst::MemFirstModelRunner bm;

	//bm.runConfig = readRunConfigFromList(runConfig);
	//bm.config = readModelConfigFromList(modelConfig);
	bm.runConfig.setFromList(runConfig);
	bm.config.setFromList(modelConfig);

	bool verbose = bm.runConfig.verbose;

	// Seed the gibbs RNG from the R RNG.
	// This allows the gibbs RNG to track with the R RNG (otherwise the gibbs random number sequence would not be reproduced).
	unsigned int uintmax = std::numeric_limits<unsigned int>::max();
	unsigned int rngSeed = uintmax * CatCont::uniformDeviate(0, 1);
	//bm.gibbs.getGenerator().seed(rngSeed);
	//bm.gibbs.iterationsPerStatusUpdate = bm.runConfig.iterationsPerStatusUpdate;
	bm.gibbs.setup(bm.runConfig.iterationsPerStatusUpdate, rngSeed);

	if (bm.config.dataType == CatCont::DataType::Circular) {
		CatCont::conditionalConfigureVMLut(bm.config.sdRanges.maxPrecision, VON_MISES_STEP_SIZE, verbose);
	}


	if (verbose) {
		Rcpp::Rcout << "Reading data." << endl;
	}
	//bm.setData(CatCont::getParticipantData(data, bm.config.dataType, verbose));
	//bm.data.setFromDataFrame(data, bm.config.dataType, verbose);
	bm.data.setup(data, bm.config.dataType, verbose);

	// Most overrides are part of model config, but equalityConstraints are passed directly.
	bm.config.overrides.equalityConstraints = convertListToMap<string>(equalityConstraints);

	// Create parameters
	if (verbose) {
		Rcpp::Rcout << "Doing parameter setup." << endl;
	}
	bm.createParameters();

	// Run the Gibbs sampler
	if (verbose) {
		Rcpp::Rcout << "Running Gibbs sampler." << endl;
	}
	bm.gibbs.run(bm.runConfig.iterations, true);


	// Collect Gibbs sampler output
	if (verbose) {
		Rcpp::Rcout << "Collecting output." << endl;
	}
	Rcpp::List priors = Rcpp::wrap(bm.priors);
	Rcpp::List posteriors = bm.gibbs.getPosteriors();


	// Save MH results
	Rcpp::List mhTuning = Rcpp::wrap(bm.mhTuningSd);
	Rcpp::DataFrame acceptanceRates = bm.gibbs.getAcceptanceRates();
	
	Rcpp::List mhList;
	mhList["tuning"] = mhTuning;
	mhList["acceptance"] = acceptanceRates;

	// Get participant numbers
	std::vector<std::string> allPnums(bm.data.participants.size());
	for (unsigned int i = 0; i < allPnums.size(); i++) {
		allPnums[i] = bm.data.participants[i].pnum;
	}

	Rcpp::List rval = Rcpp::List::create(
		Rcpp::Named("priors") = priors,
		Rcpp::Named("posteriors") = posteriors,
		Rcpp::Named("MH") = mhList,
		Rcpp::Named("pnums") = Rcpp::wrap(allPnums));

	if (verbose) {
		Rcpp::Rcout << "Done!" << endl;
	}

	return rval;
}




// [[Rcpp::export]]
Rcpp::DataFrame CCM_CPP_calculateWAIC(Rcpp::List resultsObject) {

	using namespace CatCont;

	Rcpp::List configList = resultsObject["config"];
	//Rcpp::List runConfigList = resultsObject["runConfig"];
	Rcpp::DataFrame dataDF(resultsObject["data"]);
	Rcpp::List posteriors = resultsObject["posteriors"];

	//RunConfig runCfg;
	//runCfg.setFromList(resultsObject["runConfig"]);

	ModelConfiguration config; // = readModelConfigFromList(configList);
	config.setFromList(configList);

	//vector<ParticipantData> partData = getParticipantData(data, config.dataType, false);
	DataCollection data;
	//data.setFromDataFrame(dataDF, config.dataType, false);
	data.setup(dataDF, config.dataType, false);

	if (config.dataType == CatCont::DataType::Circular) {
		conditionalConfigureVMLut(config.sdRanges.maxPrecision, VON_MISES_STEP_SIZE, false);
	}

	//CatCont::message("Copying posterior iterations");
	
	vector<ParameterList> posteriorIterations = transposePosteriors(posteriors);
	
	/*
	vector< ParameterList > posteriorIterations(runCfg.iterations);

	vector<string> posteriorNames = posteriors.names();

	for (unsigned int i = 0; i < posteriorNames.size(); i++) {

		string n = posteriorNames[i];
		Rcpp::NumericVector obs = posteriors[n];

		for (unsigned int j = 0; j < runCfg.iterations; j++) {
			posteriorIterations[j][n] = obs[j];
		}

		//CatCont::message("Done with " + n);
	}
	*/
	
	//CatCont::message("About to calculate WAIC");

	map<string, map<string, double>> waicData = CatCont::MemFirst::calculateWAIC(data, config, posteriorIterations);

	//CatCont::message("Reformatting WAIC results");

	vector<string> pnums;
	vector<double> WAIC_1_part;
	vector<double> WAIC_2_part;
	vector<double> P_1_part;
	vector<double> P_2_part;
	vector<double> LPPD_part;

	for (auto it = waicData.begin(); it != waicData.end(); it++) {

		pnums.push_back(it->first);
		map<string, double> temp = it->second;

		WAIC_1_part.push_back(temp["WAIC_1"]);
		WAIC_2_part.push_back(temp["WAIC_2"]);
		P_1_part.push_back(temp["P_1"]);
		P_2_part.push_back(temp["P_2"]);
		LPPD_part.push_back(temp["LPPD"]);

	}

	Rcpp::DataFrame result = Rcpp::DataFrame::create(Rcpp::Named("pnum") = Rcpp::wrap(pnums),
		Rcpp::Named("WAIC_1") = Rcpp::wrap(WAIC_1_part),
		Rcpp::Named("WAIC_2") = Rcpp::wrap(WAIC_2_part),
		Rcpp::Named("P_1") = Rcpp::wrap(P_1_part),
		Rcpp::Named("P_2") = Rcpp::wrap(P_2_part),
		Rcpp::Named("LPPD") = Rcpp::wrap(LPPD_part)
	);

	return result;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix CCM_CPP_MF_batchLikelihood(Rcpp::List modelConfigList, Rcpp::DataFrame dataDF, Rcpp::NumericMatrix paramMat) {
	using namespace Rcpp;
	using namespace CatCont;

	ModelConfiguration modelConfig;
	modelConfig.setFromList(modelConfigList);

	DataCollection data;
	data.setup(dataDF, modelConfig.dataType, false);

	if (modelConfig.dataType == DataType::Circular) {
		conditionalConfigureVMLut(modelConfig.sdRanges.maxPrecision, VON_MISES_STEP_SIZE, false);
	}

	CharacterVector paramNamesCV = colnames(paramMat);
	std::vector<std::string> paramNames(paramNamesCV.begin(), paramNamesCV.end());

	// 1 likelihood for each observation and each iteration.
	NumericMatrix rval(dataDF.nrow(), paramMat.nrow());

	for (size_t iter = 0; iter < (size_t)paramMat.nrow(); iter++) {
		NumericVector iterParamNV = paramMat(iter, _);
		std::vector<double> iterParamVals(iterParamNV.begin(), iterParamNV.end());
		ParameterList iterParam;
		for (size_t i = 0; i < paramNames.size(); i++) {
			iterParam[paramNames[i]] = iterParamVals[i];
		}

		//paramCont.setValues(iterParam);

		for (size_t i = 0; i < data.participants.size(); i++) {
			const string pnum = data.participants[i].pnum;

			for (size_t j = 0; j < data.conditionNames.size(); j++) {
				const ConditionData& condData = data.participants[i].condData[j];
				if (condData.study.size() == 0) {
					continue;
				}

				MemFirst::ParticipantParameters partParam = MemFirst::getParticipantParameters(iterParam, pnum, modelConfig.maxCategories);
				MemFirst::ConditionParameters condParam = MemFirst::getConditionParameters(iterParam, condData.condition);
				MemFirst::CombinedParameters combinedParam = MemFirst::combineParameters(partParam, condParam, modelConfig.sdRanges, modelConfig.dataType);

				//For each iteration, get the likelihoods for all observation for this participant/condition pair
				vector<double> likelihoods;
				if (modelConfig.dataType == DataType::Circular) {
					likelihoods = MemFirst::Circular::betweenAndWithinLikelihood(combinedParam, condData, modelConfig);
				}
				else if (modelConfig.dataType == DataType::Linear) {
					likelihoods = MemFirst::Linear::betweenAndWithinLikelihood(combinedParam, condData, modelConfig);
				}

				for (size_t obs = 0; obs < condData.originalRow.size(); obs++) {
					rval(condData.originalRow[obs], iter) = likelihoods[obs];
				}

			} // j
		} // i
	} // iter

	return rval;
}

// [[Rcpp::export]]
Rcpp::List CCM_CPP_likelihoodWrapper(Rcpp::List param, Rcpp::DataFrame data, Rcpp::List configList) {

	CatCont::ModelConfiguration modelConfig; // = readModelConfigFromList(configList);
	modelConfig.setFromList(configList);

	//CatCont::ModelVariant modelVariant = CatCont::modelVariantFromString(configList["modelVariant"]);
	//CatCont::DataType dataType = CatCont::dataTypeFromString(configList["dataType"]);


	CatCont::MemFirst::CombinedParameters cp;
	cp.pMem = param["pMem"];
	cp.pBetween = param["pBetween"];
	cp.pContBetween = param["pContBetween"];
	cp.pContWithin = param["pContWithin"];
	cp.pCatGuess = param["pCatGuess"];

	cp.contSD = param["contSD"];
	cp.cat.SD = param["catSD"];
	cp.cat.selectivity = param["catSelectivity"];

	Rcpp::NumericVector cms = param["catMu"];
	cp.cat.mu = vector<double>(cms.begin(), cms.end());


	CatCont::ConditionData cd;
	Rcpp::NumericVector study = data["study"];
	cd.study = vector<double>(study.begin(), study.end());
	Rcpp::NumericVector response = data["response"];
	cd.response = vector<double>(response.begin(), response.end());

	vector<double> likelihoods;

	if (modelConfig.dataType == CatCont::DataType::Linear) {

		//CatCont::ModelConfiguration modelConfig;
		//modelConfig.modelVariant = modelVariant;
		//modelConfig.linearConfiguration = getLinearConfigurationFromList(config);

		likelihoods = CatCont::MemFirst::Linear::betweenAndWithinLikelihood(cp, cd, modelConfig);

	} else if (modelConfig.dataType == CatCont::DataType::Circular) {

		double maxPrecision = CatCont::Circular::sdDeg_to_precRad(modelConfig.sdRanges.minSd);
		CatCont::conditionalConfigureVMLut(maxPrecision, VON_MISES_STEP_SIZE, false);

		//Parameters from degrees to radians
		cp.contSD = CatCont::Circular::sdDeg_to_precRad(cp.contSD);
		cp.cat.SD = CatCont::Circular::sdDeg_to_precRad(cp.cat.SD);
		cp.cat.selectivity = CatCont::Circular::sdDeg_to_precRad(cp.cat.selectivity);

		for (unsigned int i = 0; i < cp.cat.mu.size(); i++) {
			cp.cat.mu[i] = CatCont::Circular::degreesToRadians(cp.cat.mu[i]);
		}

		//Data from degrees to radians
		cd.study = CatCont::Circular::degreesToRadians(cd.study);
		cd.response = CatCont::Circular::degreesToRadians(cd.response);

		likelihoods = CatCont::MemFirst::Circular::betweenAndWithinLikelihood(cp, cd, modelConfig);
	}

	Rcpp::List rval = Rcpp::List::create(Rcpp::Named("likelihoods") = likelihoods);

	return rval;
}




//' Circular Distance Between Angle Vectors
//'
//' @param xs First vector of angles.
//' @param ys Second vector of angles.
//' @param absDist If `TRUE`, absolute distance is calculated. If `FALSE`, signed distance is calculated.
//' @param degrees Should be `TRUE` if angles are in degrees or `FALSE` if angles are in radians.
//' @return A vector of distances between `xs` and `ys`.
//'
//' @export
// [[Rcpp::export]]
std::vector<double> circDist(std::vector<double> xs, std::vector<double> ys, bool absDist = false, bool degrees = true) {

	if (xs.size() != ys.size()) {
		Rcpp::stop("xs and ys must have the same length.");
	}

	std::vector<double> rval(xs.size());

	for (size_t i = 0; i < rval.size(); i++) {
		rval[i] = CatCont::Circular::circularDistance(xs[i], ys[i], absDist, degrees);
	}
	
	return rval;
}

//' Clamp an Angle to a Range
//' 
//' Angles are effectively bounded in the circular space because each rotation returns to the same values.
//' Most ways of using angles, however, require treating those angles as though they were linear.
//' This function helps by restricting angles to be within a single rotation of the circle, either
//' between 0 and 360 degrees, if `pm180` is `FALSE`, or between -180 and 180 if `pm180` is `TRUE`.
//'
//' @param xs Vector of angles to clamp.
//' @param pm180 If `TRUE`, angles are clamped to the interval [-180, 180). If `FALSE`, angles are clamped to [0, 360).
//' @param degrees Should be `TRUE` if angles are in degrees or `FALSE` if angles are in radians.
//' @return A vector of clamped angles.
//'
//' @export
// [[Rcpp::export]]
std::vector<double> clampAngle(std::vector<double> xs, bool pm180 = false, bool degrees = true) {

	for (size_t i = 0; i < xs.size(); i++) {
		xs[i] = CatCont::Circular::clampAngle(xs[i], pm180, degrees);
	}

	return xs;
}


//' Get Default Priors
//' 
//' Returns default values that define the top-level prior distributions. See the "Priors" section of the package manual (Introduction.pdf).
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List getDefaultPriors() {

	map<string, double> priors = CatCont::MemFirst::MemFirstModelRunner::getDefaultPriors();
	Rcpp::List rval;
	for (auto kv : priors) {
		rval[kv.first] = kv.second;
	}

	return rval;
}

//' Get Default Metropolis-Hastings Tuning Values
//' 
//' Returns default values that affect MH update steps and the MH acceptance rate.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List getDefaultMHTuning() {

	map<string, double> mht = CatCont::MemFirst::MemFirstModelRunner::getDefaultMHTuning();
	Rcpp::List rval;
	for (auto kv : mht) {
		rval[kv.first] = kv.second;
	}

	return rval;
}


/*
The next set of functions are exported from c++ to R, but not into the package namespace for users.
*/

// This is used in an R function.
// weights should be same length as xs and sum(weights) should be 1 (enforced in R)
// [[Rcpp::export]]
double CCM_CPP_circMean(std::vector<double> xs, std::vector<double> weights, bool degrees = true) {

	if (degrees) {
		xs = CatCont::Circular::degreesToRadians(xs);
	}

	double m = CatCont::Circular::circularMean(xs, weights);

	if (degrees) {
		m = CatCont::Circular::radiansToDegrees(m);
	}

	return m;
}

// [[Rcpp::export]]
double CCM_CPP_fmod(double x, double y) {
	return std::fmod(x, y);
}


// [[Rcpp::export]]
std::vector<double> CCM_CPP_dVonMises(std::vector<double> xs, double mu, double sd, bool log, bool degrees = true, bool useLUT = true) {

	double kappa = sd;

	if (degrees) {
		xs = CatCont::Circular::degreesToRadians(xs);
		mu = CatCont::Circular::degreesToRadians(mu);
		kappa = CatCont::Circular::sdDeg_to_precRad(sd);
	}

	// Choose how to calculate dvm
	
	//CatCont::VonMisesLUT& lut = CatCont::vmLut; // Default to using vmLut
	if (useLUT) {
		if (!CatCont::vmLut.readyToUseLUT(kappa, VON_MISES_STEP_SIZE)) {
			CatCont::conditionalConfigureVMLut(kappa, VON_MISES_STEP_SIZE, true); // DEBUG
		}
	}
	else {
		CatCont::vmLut.setup(&CatCont::curriedBessel_i);
	}

	/*
	CatCont::VonMisesLUT noLut;
	if (!useLUT || !lut.ready(kappa, VON_MISES_STEP_SIZE)) {
		// Use a LUT object without the LUT, just use the bessel function.
		noLut.setup(&CatCont::curriedBessel_i);
		lut = noLut;
	}
	*/

	std::vector<double> rval(xs.size());

	for (size_t i = 0; i < xs.size(); i++) {
		double dens = CatCont::vmLut.dVonMises(xs[i], mu, kappa);

		if (degrees) {
			dens *= PI / 180.0;
		}

		if (log) {
			dens = std::log(dens);
		}
		rval[i] = dens;
	}

	return rval;
}


#ifdef USING_CAT_SPLINE

// TODO: All of the export directives here are intentionally broken.
// Replace BROKEN_EXPORT with [[Rcpp::export]]

// BROKEN_EXPORT
std::vector<double> dZeroDerivSpline(std::vector<double> zs) {
	std::vector<double> rval(zs.size());

	for (size_t i = 0; i < zs.size(); i++) {
		rval[i] = CatCont::zeroDerivativeCubicSplineDensity(zs[i]);
	}

	return rval;
}

// BROKEN_EXPORT
double dPlatSplineFull_R(double x, double mu, double platHW, double splineHW, bool linear) {
	return CatCont::dPlatSplineFull(x, mu, platHW, splineHW, linear, true); // degrees = true
}

// BROKEN_EXPORT
double dPlatSpline_R(double absDist, double platHW, double splineHW) {
	return CatCont::dPlatSpline(absDist, platHW, splineHW);
}

// BROKEN_EXPORT
std::vector<double> dPlatSplineWeights_R(double x, std::vector<double> mu, std::vector<double> platHW, double splineHW, bool linear = false) {
	CatCont::CategoryParameters catPar;
	catPar.mu = mu;
	//catPar.SD = catSD;
	catPar.platHW = platHW;
	catPar.selectivity = splineHW;

	CatCont::ModelConfiguration modCfg;
	modCfg.dataType = linear ? CatCont::DataType::Linear : CatCont::DataType::Circular;
	
	return CatCont::dPlatSplineWeights(x, catPar, modCfg);
}

// BROKEN_EXPORT
std::vector<double> dPlatSpline_WC(double x, std::vector<double> mu, std::vector<double> platHW, double splineHW, bool linear = false) {

	// TODO: Make splineHW a vector

	if (mu.size() != platHW.size()) {
		Rcpp::stop("length(mu) != length(platHW).");
	}

	CatCont::ModelConfiguration modCfg;
	modCfg.weightsDistribution = CatCont::WeightsDistribution::PlatSpline;
	modCfg.dataType = linear ? CatCont::DataType::Linear : CatCont::DataType::Circular;
	modCfg.maxCategories = mu.size();

	CatCont::CategoryParameters catPar;
	if (linear) {
		catPar.mu = mu;
		//catPar.SD = catSD;
		catPar.platHW = platHW;
		catPar.selectivity = splineHW;
	} else {
		// circular
		x = CatCont::Circular::degreesToRadians(x);

		catPar.mu = CatCont::Circular::degreesToRadians(mu);
		//catPar.SD = catSD;
		catPar.platHW = CatCont::Circular::degreesToRadians(platHW);

		// If using PlatSpline, you don't want precision, but do want to convert from deg to rad.
		catPar.selectivity = CatCont::Circular::degreesToRadians(splineHW);
		//catPar.selectivity = CatCont::Circular::sdDeg_to_precRad(splineHW);
	}

	CatCont::WeightsCalculator wc(&modCfg);

	wc.calcWeights(x, catPar);

	// This is ok because weights.size() == modCfg.maxCategories == mu.size() for this function
	return wc.weights; 
}

// BROKEN_EXPORT
double calcPlatSplineLambda_R(std::vector<double> weights, std::string lambdaVariant) {

	CatCont::ModelConfiguration modCfg;
	modCfg.weightsDistribution = CatCont::WeightsDistribution::PlatSpline;
	//modCfg.dataType = linear ? CatCont::DataType::Linear : CatCont::DataType::Circular; // doesn't matter for lambda
	modCfg.lambdaVariant = CatCont::lambdaVariantFromString(lambdaVariant);

	CatCont::WeightsCalculator wc(&modCfg);
	wc.weights = weights;
	wc.weightCount = weights.size();

	//return CatCont::calcPlatSplineLambda(weights, CatCont::LambdaVariant::CatWeightSum);

	return wc.calcLambda();
}

#endif // USING_CAT_SPLINE


/*
namespace CatCont {

	vector<ParticipantData> getParticipantData(Rcpp::DataFrame df, CatCont::DataType dataType, bool verbose) {

		Rcpp::CharacterVector pnumsColRaw = df["pnum"];
		vector<string> pnumsCol(pnumsColRaw.begin(), pnumsColRaw.end());

		Rcpp::CharacterVector condsColRaw = df["cond"];
		vector<string> condsCol(condsColRaw.begin(), condsColRaw.end());

		Rcpp::NumericVector studyColRaw = df["study"];
		vector<double> studyCol(studyColRaw.begin(), studyColRaw.end());

		Rcpp::NumericVector respColRaw = df["response"];
		vector<double> respCol(respColRaw.begin(), respColRaw.end());

		return copyParticipantData(pnumsCol, condsCol, studyCol, respCol, dataType, verbose);

	}

} // namespace CatCont


CatCont::Linear::LinearConfiguration getLinearConfigurationFromList(Rcpp::List configList) {

	CatCont::Linear::LinearConfiguration lc;

	Rcpp::NumericVector rr = configList["responseRange"];
	Rcpp::NumericVector cmr = configList["catMuRange"];

	lc.response.lower = rr[0];
	lc.response.upper = rr[1];

	lc.catMu.lower = cmr[0];
	lc.catMu.upper = cmr[1];

	return lc;
}


CatCont::RunConfig readRunConfigFromList(Rcpp::List runCfgList) {

	CatCont::RunConfig rval;

	rval.iterations = runCfgList["iterations"];
	rval.iterationsPerStatusUpdate = runCfgList["iterationsPerStatusUpdate"];


	if (runCfgList.containsElementNamed("verbose")) {
		int verbose_int = runCfgList["verbose"];
		rval.verbose = (verbose_int == 1);
	}

	if (runCfgList.containsElementNamed("profileParameterTypes")) {
		int profile_int = runCfgList["profileParameterTypes"];
		rval.profileParameterTypes = (profile_int == 1);
	}

	return rval;
}

CatCont::ModelConfiguration readModelConfigFromList(Rcpp::List configList) {

	CatCont::ModelConfiguration config;

	string modelVariantStr = configList["modelVariant"];
	config.modelVariant = CatCont::modelVariantFromString(modelVariantStr);

	string dataTypeStr = configList["dataType"];
	config.dataType = CatCont::dataTypeFromString(dataTypeStr);

	string ccName = configList["cornerstoneConditionName"];
	config.cornerstoneConditionName = ccName;

	config.catMuPriorApproximationPrecision = configList["catMuPriorApproximationPrecision"];

	if (config.modelVariant == CatCont::ModelVariant::ZL) {
		config.maxCategories = 0;
	}
	else {
		config.maxCategories = configList["maxCategories"];
	}

	if (config.dataType == CatCont::DataType::Linear) {
		config.linearConfiguration = getLinearConfigurationFromList(configList);
	}

	if (configList.containsElementNamed("calculateParticipantLikelihoods")) {
		// Going straight to bool doesn't work right
		int calculateParticipantLikelihoods_temp = configList["calculateParticipantLikelihoods"];
		config.calculateParticipantLikelihoods = (calculateParticipantLikelihoods_temp == 1);
	}

	// Set ranges
	config.sdRanges.minSd = configList["minSD"];
	config.sdRanges.maxSd = VON_MISES_MAX_SD;
	config.sdRanges.minPrecision = CatCont::Circular::sdDeg_to_precRad(config.sdRanges.maxSd);
	config.sdRanges.maxPrecision = CatCont::Circular::sdDeg_to_precRad(config.sdRanges.minSd);


	// Condition Effects
	Rcpp::List conditionEffects = configList["conditionEffects"];
	Rcpp::CharacterVector rawNames = conditionEffects.names();
	for (int i = 0; i < rawNames.size(); i++) {
		string name = Rcpp::as<string>(rawNames[i]);
		Rcpp::CharacterVector ce = conditionEffects[name];

		vector<string> cev(ce.size());

		for (int j = 0; j < ce.size(); j++) {
			cev[j] = Rcpp::as<string>(ce[j]);
		}

		config.conditionEffects[name] = cev;
	}

	for (auto it = config.conditionEffects.begin(); it != config.conditionEffects.end(); it++) {
		const vector<string>& ce = it->second;

		bool isNone = ce.size() == 0 || (ce.size() == 1 && ce[0] == "none");

		if (!isNone) {
			config.paramWithConditionEffects.push_back(it->first);
		}
	}

	
	//Rcpp::CharacterVector paramWithConditionEffects = configList["parametersWithConditionEffects"];
	//for (unsigned int i = 0; i < (unsigned int)paramWithConditionEffects.size(); i++) {
	//	config.paramWithConditionEffects.push_back((string)paramWithConditionEffects[i]);
	//}
	

	// Overrides are now part of the config list
	if (configList.containsElementNamed("priorOverrides")) {
		config.overrides.priors = convertListToMap<double>(configList["priorOverrides"]);
	}
	if (configList.containsElementNamed("mhTuningOverrides")) {
		config.overrides.mhTunings = convertListToMap<double>(configList["mhTuningOverrides"]);
	}
	if (configList.containsElementNamed("startingParamValues")) {
		config.overrides.startingValues = convertListToMap<double>(configList["startingParamValues"]);
	}
	if (configList.containsElementNamed("constantParamValues")) {
		config.overrides.constantValues = convertListToMap<double>(configList["constantParamValues"]);
	}

	//if (configList.containsElementNamed("equalityConstraints")) {
	//	config.overrides.equalityConstraints = convertListToMap<double>(configList["equalityConstraints"]);
	//}


	// privateConfig
	if (configList.containsElementNamed("privateConfig")) {
		Rcpp::List pcList = configList["privateConfig"];

		// TODO: Where does this go? It doesn't appear to be used anywhere (why not?)
		if (pcList.containsElementNamed("useVonMisesLookupTable")) {
			int vmlut_int = pcList["useVonMisesLookupTable"];
			config.privateConfig.useVonMisesLookupTable = (vmlut_int == 1);
		}

		// Temporarily, these settings are passed to c++ in privateConfig
		if (pcList.containsElementNamed("catMuShared")) {
			int cms_int = pcList["catMuShared"];
			config.catMuShared = (cms_int == 1);
		}
		if (pcList.containsElementNamed("catActiveShared")) {
			int cas_int = pcList["catActiveShared"];
			config.catActiveShared = (cas_int == 1);
		}
		if (pcList.containsElementNamed("catActiveHierPrior")) {
			int cahp_int = pcList["catActiveHierPrior"];
			config.catActiveHierPrior = (cahp_int == 1);
		}
		if (pcList.containsElementNamed("catActiveDistancePrior")) {
			int cadp_int = pcList["catActiveDistancePrior"];
			config.catActiveDistancePrior = (cadp_int == 1);
		}

		// If catMu not shared, catActive can't be shared
		if (!config.catMuShared) {
			config.catActiveShared = false;
		}

		// Can't have hier prior if shared
		if (config.catActiveShared) {
			config.catActiveHierPrior = false;
		}

		//
		if (pcList.containsElementNamed("withinItem_contSD_is_withinSD")) {
			int bool_int = pcList["withinItem_contSD_is_withinSD"];
			config.withinItem_contSD_is_withinSD = (bool_int == 1);
		}
	}

#ifdef USING_CAT_SPLINE
	if (configList.containsElementNamed("weightsDistribution")) {
		string weightsDistributionStr = configList["weightsDistribution"];
		config.weightsDistribution = CatCont::weightsDistributionFromString(weightsDistributionStr);
	}
	else {
		config.weightsDistribution = CatCont::WeightsDistribution::Default;
	}
	config.lambdaVariant = CatCont::LambdaVariant::None; // TODO: Choose value.
#endif

	return config;
}
*/


#endif // COMPILING_WITH_RCPP