#include "CCM_BayesianModel.h"


#ifdef COMPILING_WITH_RCPP

namespace CatCont {

void Bayesian::_doMhOverrides(void) {
	Rcpp::Rcout << "MH overrides: " << rcppConfig.mhTuningOverrides.size() << endl;
	if (rcppConfig.mhTuningOverrides.size() == 0) {
		return;
	}

	Rcpp::CharacterVector rawNames = rcppConfig.mhTuningOverrides.names();
	std::vector<string> names(rawNames.begin(), rawNames.end());

	for (unsigned int i = 0; i < names.size(); i++) {

		string name = names[i];

		if (mhTuningSd.find(name) != mhTuningSd.end()) {

			Rcpp::NumericVector vVal = rcppConfig.mhTuningOverrides[name];
			double val = vVal[0];

			mhTuningSd[name] = val;
		} else {
			logMessage("_doMhOverrides", "Warning: Invalid MH tuning override key: \"" + names[i] + "\" ignored.");
		}

	}
}

void Bayesian::_doPriorOverrides(void) {
	Rcpp::Rcout << "Prior overrides: " << rcppConfig.priorOverrides.size() << endl;
	if (rcppConfig.priorOverrides.size() == 0) {
		return;
	}

	Rcpp::CharacterVector rawNames = rcppConfig.priorOverrides.names();
	std::vector<string> names(rawNames.begin(), rawNames.end());

	for (unsigned int i = 0; i < names.size(); i++) {

		string name = names[i];

		if (priors.find(name) != priors.end()) {

			Rcpp::NumericVector vVal = rcppConfig.priorOverrides[name];
			double val = vVal[0];

			priors[name] = val;
		} else {
			logMessage("_doPriorOverrides", "Warning: Invalid prior override key: \"" + name + "\" ignored.");
		}

	}

}

void Bayesian::_doConstantParameterOverrides(void) {

	Rcpp::Rcout << "Constant parameter value overrides: " << rcppConfig.constantValueOverrides.size() << endl;
	if (rcppConfig.constantValueOverrides.size() == 0) {
		return;
	}

	Rcpp::CharacterVector rawNames = rcppConfig.constantValueOverrides.names();
	std::vector<string> names(rawNames.begin(), rawNames.end());
	map<string, double> vals;
	for (unsigned int i = 0; i < names.size(); i++) {
		vals[names[i]] = rcppConfig.constantValueOverrides[names[i]];
	}

	vector<GibbsParameter*> parameters = gibbs.getParameters();
	for (unsigned int i = 0; i < parameters.size(); i++) {
		GibbsParameter* par = parameters[i];

		map<string, double>::iterator it = vals.find(par->name);

		if (it != vals.end()) {
			gibbs.replaceParameter(par->name, ConstantParameter(it->second, par->name, par->group));
		}

	}

	//Also check which values are provided but for which there is no parameter.
	for (map<string, double>::iterator it = vals.begin(); it != vals.end(); it++) {
		if (!gibbs.hasParameter(it->first)) {
			std::stringstream ss;
			ss << "Note: Constant value provided for \"" << it->first << "\", but there is no parameter by that name.";
			logMessage("setParameterStartingValues", ss.str());
		}
	}

}

void Bayesian::_doStartingValueOverrides(void) {
	Rcpp::Rcout << "Starting value overrides: " << rcppConfig.startingValueOverrides.size() << endl;
	if (rcppConfig.startingValueOverrides.size() == 0) {
		return;
	}

	Rcpp::CharacterVector rawNames = rcppConfig.startingValueOverrides.names();
	std::vector<string> names(rawNames.begin(), rawNames.end());

	std::map<std::string, double> startingValues;
	for (unsigned int i = 0; i < names.size(); i++) {
		startingValues[names[i]] = rcppConfig.startingValueOverrides[names[i]];
	}

	this->setParameterStartingValues(startingValues);
}


vector<ParticipantData> getParticipantData(Rcpp::DataFrame df, CatCont::DataType dataType, bool verbose) {

	Rcpp::CharacterVector pnumsColRaw = df["pnum"];
	vector<string> pnumsCol(pnumsColRaw.begin(), pnumsColRaw.end());

	Rcpp::CharacterVector condsColRaw = df["cond"];
	vector<string> condsCol(condsColRaw.begin(), condsColRaw.end());

	Rcpp::NumericVector studyColRaw = df["study"];
	vector<double> studyCol(studyColRaw.begin(), studyColRaw.end());

	Rcpp::NumericVector respColRaw = df["response"];
	vector<double> respCol(respColRaw.begin(), respColRaw.end());


	Rcpp::CharacterVector uniquePnums_raw = Rcpp::unique(pnumsColRaw).sort();
	vector<string> uniquePnums(uniquePnums_raw.begin(), uniquePnums_raw.end());

	Rcpp::CharacterVector uniqueConditions_raw = Rcpp::unique(condsColRaw).sort();
	vector<string> uniqueConditions(uniqueConditions_raw.begin(), uniqueConditions_raw.end());

	vector<ParticipantData> data;

	for (unsigned int p = 0; p < uniquePnums.size(); p++) {

		ParticipantData thisPart;
		thisPart.pnum = uniquePnums[p];

		//Find which rows of the data correspond to the given pnum
		vector<unsigned int> pnumRows;
		for (unsigned int row = 0; row < pnumsCol.size(); row++) {
			if (pnumsCol[row] == uniquePnums[p]) {
				pnumRows.push_back(row);
			}
		}

		for (unsigned int c = 0; c < uniqueConditions.size(); c++) {

			//Copy from the secondary data frame to primitive data types in a ConditionData struct
			ConditionData condData;

			condData.condition = uniqueConditions[c];

			for (unsigned int i = 0; i < pnumRows.size(); i++) {

				unsigned int thisRow = pnumRows[i];

				//If the condition on the current row is equal to the cth condition, store the data
				if (condsCol[thisRow] == uniqueConditions[c]) {
					double study = studyCol[thisRow];
					double response = respCol[thisRow];

					if (dataType == DataType::Circular) {
						study = Circular::degreesToRadians(study); //TODO: This should probably be done by the model...
						response = Circular::degreesToRadians(response);
					}
					condData.study.push_back(study);
					condData.response.push_back(response);
				}
			}

			thisPart.condData.push_back(condData);

			if (verbose) {
				Rcpp::Rcout << "Participant " << uniquePnums[p] << ", condition " << uniqueConditions[c] << ": " <<
					_convertToString(condData.study.size()) << " trials found." << endl;
			}
		}

		data.push_back(thisPart);

	}

	return data;

}

} // namespace CatCont

CatCont::Linear::LinearConfiguration getLinearConfigurationFromList(Rcpp::List configList) {

	CatCont::Linear::LinearConfiguration lc;

	string modelVariantStr = configList["modelVariant"];
	lc.modelVariant = CatCont::modelVariantFromString(modelVariantStr);

	Rcpp::NumericVector rr = configList["responseRange"];
	Rcpp::NumericVector cmr = configList["catMuRange"];

	lc.response.lower = rr[0];
	lc.response.upper = rr[1];

	lc.catMu.lower = cmr[0];
	lc.catMu.upper = cmr[1];

	return lc;
}

CatCont::Bayesian::Configuration readConfigurationFromList(Rcpp::List configList) {

	CatCont::Bayesian::Configuration config;

	config.iterations = configList["iterations"];
	config.maxCategories = configList["maxCategories"];

	string ccName = configList["cornerstoneConditionName"];
	config.cornerstoneConditionName = ccName;

	string modelVariantStr = configList["modelVariant"];
	config.modelVariant = CatCont::modelVariantFromString(modelVariantStr);

	if (config.modelVariant == CatCont::ModelVariant::ZL) {
		config.maxCategories = 0;
	}

	string dataTypeStr = configList["dataType"];
	config.dataType = CatCont::dataTypeFromString(dataTypeStr);

	if (config.dataType == CatCont::DataType::Linear) {
		config.linearConfiguration = getLinearConfigurationFromList(configList);
	}


	//Going straight to bool doesn't work right
	int calculateParticipantLikelihoods_temp = configList["calculateParticipantLikelihoods"];
	config.calculateParticipantLikelihoods = (calculateParticipantLikelihoods_temp == 1);

	config.ranges.minSd = configList["minSD"];
	config.ranges.maxSd = VON_MISES_MAX_SD;
	config.ranges.minPrecision = CatCont::Circular::sdDeg_to_precRad(config.ranges.maxSd);
	config.ranges.maxPrecision = CatCont::Circular::sdDeg_to_precRad(config.ranges.minSd);


	Rcpp::CharacterVector paramWithConditionEffects = configList["parametersWithConditionEffects"];
	for (unsigned int i = 0; i < (unsigned int)paramWithConditionEffects.size(); i++) {
		config.paramWithConditionEffects.push_back((string)paramWithConditionEffects[i]);
	}

	return config;
}

double curriedBesselFunction(double x) {
	return R::bessel_i(x, 0, 2); //Exponent scale == true (2 == true, lol)
}

void conditionalConfigureVMLut(double maxValue, double stepSize) {

	bool rangeCorrect = abs(CatCont::vmLut.maxValue() - maxValue) < 0.00001;
	bool stepSizeCorrect = abs(CatCont::vmLut.stepSize - stepSize) < 0.00001;

	if (rangeCorrect && stepSizeCorrect) {
		Rcpp::Rcout << "Von Mises look up table already set up." << endl;
	} else {
		Rcpp::Rcout << "Setting up Von Mises look up table." << endl;
		CatCont::vmLut.setup(maxValue, stepSize, &curriedBesselFunction);
	}
}

// [[Rcpp::export]]
Rcpp::DataFrame CCM_CPP_calculateWAIC(Rcpp::List resultsObject) {

	using namespace CatCont;

	Rcpp::List configList = resultsObject["config"];
	Rcpp::DataFrame data = resultsObject["data"];
	Rcpp::List posteriors = resultsObject["posteriors"];

	Bayesian::Configuration config = readConfigurationFromList(configList);

	vector<ParticipantData> partData = getParticipantData(data, config.dataType, false);

	if (config.dataType == CatCont::DataType::Circular) {
		conditionalConfigureVMLut(config.ranges.maxPrecision, VON_MISES_STEP_SIZE);
	}

	vector< ParameterList > posteriorIterations(config.iterations);

	vector<string> posteriorNames = posteriors.names();

	for (unsigned int i = 0; i < posteriorNames.size(); i++) {

		string n = posteriorNames[i];
		Rcpp::NumericVector obs = posteriors[n];

		for (unsigned int j = 0; j < config.iterations; j++) {
			posteriorIterations[j][n] = obs[j];
		}
	}

	Rcpp::Rcout << "Calculating WAIC" << std::endl;

	map<string, map<string, double>> waicData = CatCont::calculateWAIC(partData, config, posteriorIterations);

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

	//return Rcpp::wrap(waicData);
}




// [[Rcpp::export]]
Rcpp::List CCM_CPP_runParameterEstimation(Rcpp::List generalConfig,
	Rcpp::DataFrame data,
	Rcpp::List mhTuningOverrides,
	Rcpp::List priorOverrides,
	Rcpp::List startingValueOverrides,
	Rcpp::List constantValueOverrides
)
{
	CatCont::Bayesian bm;

	bm.config = readConfigurationFromList(generalConfig);

	bm.gibbs.iterationsPerStatusUpdate = generalConfig["iterationsPerStatusUpdate"];

	if (bm.config.dataType == CatCont::DataType::Circular) {
		conditionalConfigureVMLut(bm.config.ranges.maxPrecision, VON_MISES_STEP_SIZE);
	}


	Rcpp::Rcout << "Reading data." << endl;

	bm.setData(CatCont::getParticipantData(data, bm.config.dataType, true));

	bm.rcppConfig.mhTuningOverrides = mhTuningOverrides;
	bm.rcppConfig.priorOverrides = priorOverrides;
	bm.rcppConfig.startingValueOverrides = startingValueOverrides;
	bm.rcppConfig.constantValueOverrides = constantValueOverrides;


	Rcpp::Rcout << "Doing parameter setup." << endl;

	bm.setMhTuning();
	bm.setPriors();
	bm.createParameters(); //must happen after priors and mhTuning

	Rcpp::Rcout << "Running Gibbs sampler." << endl;

	bm.gibbs.run(bm.config.iterations, true);


	Rcpp::Rcout << "Collecting output." << endl;

	//Get output from the sampler
	Rcpp::List posteriors = bm.gibbs.getPosteriors();
	Rcpp::DataFrame acceptanceRates = bm.gibbs.getAcceptanceRates();

	Rcpp::List mhTuning = Rcpp::wrap(bm.mhTuningSd);
	Rcpp::List priors = Rcpp::wrap(bm.priors);

	//Get participant numbers
	std::vector<std::string> allPnums(bm.data.participants.size());
	for (unsigned int i = 0; i < allPnums.size(); i++) {
		allPnums[i] = bm.data.participants[i].pnum;
	}

	Rcpp::List rval = Rcpp::List::create(Rcpp::Named("posteriors") = posteriors,
		Rcpp::Named("mhAcceptance") = acceptanceRates,
		Rcpp::Named("mhTuning") = mhTuning,
		Rcpp::Named("priors") = priors,
		Rcpp::Named("constantValueOverrides") = constantValueOverrides,
		Rcpp::Named("pnums") = Rcpp::wrap(allPnums),
		Rcpp::Named("conditionNames") = bm.data.conditionNames);

	return rval;
}



// [[Rcpp::export]]
Rcpp::List CCM_CPP_likelihoodWrapper(Rcpp::List param, Rcpp::DataFrame data, Rcpp::List config) {

	CatCont::ModelVariant modelVariant = CatCont::modelVariantFromString(config["modelVariant"]);
	CatCont::DataType dataType = CatCont::dataTypeFromString(config["dataType"]);


	CatCont::CombinedParameters cp;
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

	if (dataType == CatCont::DataType::Linear) {

		CatCont::Linear::LinearConfiguration lc = getLinearConfigurationFromList(config);

		likelihoods = CatCont::Linear::betweenAndWithinLikelihood(cp, cd, lc);

	} else if (dataType == CatCont::DataType::Circular) {

		double maxPrecision = CatCont::Circular::sdDeg_to_precRad(config["minSD"]);
		conditionalConfigureVMLut(maxPrecision, VON_MISES_STEP_SIZE);

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

		likelihoods = CatCont::Circular::betweenAndWithinLikelihood(cp, cd, modelVariant);
	}

	Rcpp::List rval = Rcpp::List::create(Rcpp::Named("likelihoods") = likelihoods);

	return rval;
}

#endif //COMPILING_WITH_RCPP