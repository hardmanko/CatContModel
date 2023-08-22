#include "CF_Model.h"

#include "CCM_EqualityConstraints.h"

// [[Rcpp::export]]
Rcpp::List CCM_CF_RPE(Rcpp::DataFrame dataDF, Rcpp::List modCfgList, Rcpp::List runConfigList) {
	CatCont::CatFirst::CatFirstModel mod;
	if (!mod.setup(dataDF, modCfgList, runConfigList)) {
		CatCont::stop("Model setup failed.");
	}


	// Run the Gibbs sampler
	if (mod.runConfig.verbose) {
		Rcpp::Rcout << "Running Gibbs sampler." << endl;
	}
	mod.gibbs.run(mod.runConfig.iterations, true);


	// Collect Gibbs sampler output
	if (mod.runConfig.verbose) {
		Rcpp::Rcout << "Collecting output." << endl;
	}
	Rcpp::List priors = Rcpp::wrap(mod.priors);
	Rcpp::List posteriors = mod.gibbs.getPosteriors();


	// Save MH results
	Rcpp::List mhTuning = Rcpp::wrap(mod.mhTuningSd);
	Rcpp::DataFrame acceptanceRates = mod.gibbs.getAcceptanceRates();

	Rcpp::List mhList;
	mhList["tuning"] = mhTuning;
	mhList["acceptance"] = acceptanceRates;

	// Get participant numbers in order
	std::vector<std::string> allPnums(mod.data.participants.size());
	for (size_t i = 0; i < allPnums.size(); i++) {
		allPnums[i] = mod.data.participants[i].pnum;
	}

	Rcpp::List rval = Rcpp::List::create(
		Rcpp::Named("priors") = priors,
		Rcpp::Named("posteriors") = posteriors,
		Rcpp::Named("MH") = mhList,
		Rcpp::Named("pnums") = Rcpp::wrap(allPnums));

	if (mod.runConfig.verbose) {
		Rcpp::Rcout << "Done!" << endl;
	}

	return rval;
}

namespace CatCont {
namespace CatFirst {



bool CatFirstModel::setup(Rcpp::DataFrame dataDF, Rcpp::List modCfgList, Rcpp::List runConfigList) {
	// Model config
	if (!this->config.setFromList(modCfgList)) {
		return false;
	}

	/*
	this->modelConfig.setFromList(modCfgList);

	Rcpp::List privateCfg = modCfgList["privateConfig"];

	// CatFirst specific model config
	Rcpp::List cfmcl = privateCfg["CatFirst"];
	this->CFConfig.setFromList(cfmcl);

	// Variable Precision config
	if (privateCfg.containsElementNamed("VPConfig")) {
		this->VPConfig.setFromList(privateCfg["VPConfig"]);
	}
	*/


	// Run config
	//this->runConfig = readRunConfigFromList(runConfig);
	this->runConfig.setFromList(runConfigList);

	// Set the data
	//this->data.setFromDataFrame(dataDF, config.modelConfig.dataType, true);
	this->data.setup(dataDF, config.modelConfig.dataType, true);

	// Set up likelihood calculator (which also sets up von Mises LUT)
	_likeCalc.setup(config);

	// Set priors and MH tuning
	_setPriors();
	_setMHTuning();

	// Configure the Gibbs sampler

	// Seed the gibbs RNG from the R RNG.
	// This allows the gibbs RNG to track with the R RNG (otherwise the gibbs random number sequence would not be reproduced).
	unsigned int uintmax = std::numeric_limits<unsigned int>::max();
	unsigned int rngSeed = uintmax * CatCont::uniformDeviate(0, 1);

	//this->gibbs.getGenerator().seed(rngSeed);
	//this->gibbs.iterationsPerStatusUpdate = this->runConfig.iterationsPerStatusUpdate;

	this->gibbs.setup(this->runConfig.iterationsPerStatusUpdate, rngSeed);

	if (runConfig.verbose) {
		CatCont::message("Setting up Gibbs parameters");
	}
	this->createGibbsParameters();

	if (!this->gibbs.prepareToRun()) {
		CatCont::stop("Gibbs sampler not able to run.");
	}

	// After creating parameters, setup the parameter extractor
	if (runConfig.verbose) {
		CatCont::message("Setting up parameter extractor");
	}
	ParamContainer startingVals = this->gibbs.getCurrentParameterValues();

	// DEBUG
	/*
	vector<string> pNames = startingVals.reconstructOrderedNames();
	CatCont::message("All parameter names:");
	for (const string& pn : pNames) {
		CatCont::message(pn);
	}
	*/

	_paramExtractor.setup(config, this->data, startingVals);

	// TODO: More setup?

	return true;
}


void CatFirstModel::createGibbsParameters(void) {

	using namespace std::placeholders;

	// TODO: Set priors and MH tuning. Should be in setup?

	gibbs.sectionTracker.sectionStart("Total");


	////////////////////////////////////////////////////////////////////////////
	// Combined participant and population

	gibbs.sectionTracker.sectionStart("Participant Standard");
	if (runConfig.verbose) {
		CatCont::message("-- Param: Participant Standard");
	}

	vector<string> standardParam = { "pMem", "pBetween", "pContBetween", "pContWithin", "pCatGuess", "contSD", "catSD", "catSelectivity" };
	vector<bool> isProbParam = { true, true, true, true, true, false, false, false };

	for (size_t pi = 0; pi < standardParam.size(); pi++) {

		const string& baseParamName = standardParam[pi];
		const string paramName_part = baseParamName + "_part";
		bool isShared = config.modelConfig.sharedParameters.count(baseParamName) == 1;
		//bool isShared = config.cfConfig.sharedParticipantParam.count(baseParamName) == 1;

		// Population param
		if (isShared) {
			// Mu is estimated based on all data, using .mu.mu and .mu.var as priors.
			MH_Parameter partMuShared;
			partMuShared.name = paramName_part + ".mu";
			partMuShared.group = "populationParam";

			partMuShared.llFunction = std::bind(&CatFirstModel::_LLS_PR_participantStandard_shared, this, _1, _2, priors.at(paramName_part + ".mu.mu"), priors.at(paramName_part + ".mu.var"));
			partMuShared.deviateFunction = std::bind(CatCont::normalDeviate, _1, mhTuningSd.at(paramName_part));

			double startValue = CatCont::uniformDeviate(config.modelConfig.sdRanges.minSd + 2, 40);
			gibbs.addParameter(partMuShared, startValue);

			// Variance is constant whatever (not used)
			ConstantParameter partVarShared(1, paramName_part + ".var", "populationParam");
			gibbs.addParameter(partVarShared);
		}
		else {
			ConjugateParameter popMu;

			popMu.name = paramName_part + ".mu";
			popMu.group = "populationParam";

			popMu.samplingFunction = std::bind(&CatFirstModel::_conjugateSample_mu, this, _1, baseParamName, priors.at(paramName_part + ".mu.mu"), priors.at(paramName_part + ".mu.var"));

			gibbs.addParameter(popMu, 0);


			ConjugateParameter popVar;

			popVar.name = paramName_part + ".var";
			popVar.group = "populationParam";

			popVar.samplingFunction = std::bind(&CatFirstModel::_conjugateSample_var, this, _1, baseParamName, priors.at(paramName_part + ".var.a"), priors.at(paramName_part + ".var.b"));

			gibbs.addParameter(popVar, 30);
		}

		// Participant param
		for (size_t i = 0; i < data.participants.size(); i++) {

			string indexPartParamName = paramName_part + "[" + data.participants[i].pnum + "]";

			if (isShared) {
				DependentParameter dp(indexPartParamName, paramName_part, paramName_part + ".mu");
				gibbs.addParameter(dp);
			}
			else {
				// Not shared: Estimate
				MH_Parameter par;
				par.name = indexPartParamName;
				par.group = paramName_part;

				par.llFunction = std::bind(&CatFirstModel::_LLS_PR_participantStandard, this, _1, _2, i, baseParamName);
				par.deviateFunction = std::bind(CatCont::normalDeviate, _1, mhTuningSd.at(paramName_part));

				double startValue = 1;
				if (isProbParam[pi]) {
					startValue = CatCont::uniformDeviate(-2, 2);
				}
				else {
					// assume SD param
					startValue = CatCont::uniformDeviate(config.modelConfig.sdRanges.minSd + 2, 40);
				}

				gibbs.addParameter(par, startValue);
			}
			
		} // i
	} // standardParam

	gibbs.sectionTracker.sectionEnd("Participant Standard");

	// End Combine participant with population
	////////////////////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////////////////////////
	// Participant Category
	gibbs.sectionTracker.sectionStart("Participant Category");
	if (runConfig.verbose) {
		CatCont::message("-- Param: Participant Category");
	}

	// The start value for catMu are a grid within the response range.
	// Note that this uses the data response range rather than the user-settable config response range.
	// This works in the same way for both linear and circular.
	std::vector<double> catMuStartGrid(config.modelConfig.maxCategories);
	double catMuGridStepSize = (data.responseRange.upper - data.responseRange.lower) / (double)config.modelConfig.maxCategories;
	for (size_t k = 0; k < config.modelConfig.maxCategories; k++) {
		catMuStartGrid[k] = (k + 0.5) * catMuGridStepSize + data.responseRange.lower;
	}
	// Start grid uses data range, but data are in radians and catMu param are in degrees.
	if (config.modelConfig.dataType == DataType::Circular) {
		catMuStartGrid = CatCont::Circular::radiansToDegrees(catMuStartGrid);
	}

	//const double catMuCandidateSd = mhTuningSd.at("catMu_part"); //linear
	//if (config.modelConfig.dataType == DataType::Circular) {
	//	catMuCandidateSd = CatCont::Circular::degreesToRadians(catMuCandidateSd);
	//}
	const bool shouldModuloCatMu = config.modelConfig.dataType == CatCont::DataType::Circular;

	// If a category parameter is shared between participants, 
	// use the first participant as the holder for the shared parameters.
	// Copy the shared values to all other participants with DependentParameter.
	unsigned int sharedHolderIndex = 0;

	// or this->_catMuShared etc.
	//config.modelConfig.catMuShared = config.cfConfig.sharedParticipantParam.count("catMu") == 1;
	//config.modelConfig.catActiveShared = config.cfConfig.sharedParticipantParam.count("catActive") == 1;

	bool catMuIsShared = config.modelConfig.sharedParameters.count("catMu") == 1;
	this->_catActiveShared = config.modelConfig.sharedParameters.count("catActive") == 1;
	this->_catTypeShared = config.modelConfig.sharedParameters.count("catType") == 1;
	this->_catTypeNames = config.cfConfig.catTypes.at("type");

	size_t catTypeCSIndex = config.cfConfig.getCatTypeCornerstoneIndex();

	string catMuGroup = "catMu_part"; // catMu_ik?
	string catActiveGroup = "catActive_part";

	for (size_t i = 0; i < data.participants.size(); i++) {
		for (size_t k = 0; k < config.modelConfig.maxCategories; k++) {

			string catIstr = "[" + data.participants.at(i).pnum + "," + CatCont::catIndexString(k) + "]";


			// ----------------------------------------------
			// catMu
			if (catMuIsShared) {

				if (i == sharedHolderIndex) {
					MH_Parameter sharedCatMu;
					sharedCatMu.name = "catMu_part" + catIstr;
					sharedCatMu.group = catMuGroup;

					//sharedCatMu.deviateFunction = bind(CatCont::normalDeviate, _1, catMuCandidateSd);
					sharedCatMu.deviateFunction = bind(CatCont::MemFirst::catMuDeviateFunction, _1, mhTuningSd.at("catMu_part"), shouldModuloCatMu);
					sharedCatMu.llFunction = bind(&CatFirstModel::_LLS_PR_catMu_shared, this, _1, _2, k);

					gibbs.addParameter(sharedCatMu, catMuStartGrid[k]);
				}
				else {
					// If not the holder, catMu and catActive are dependent parameters
					string catMuSource = "catMu_part[" + data.participants.at(sharedHolderIndex).pnum + "," + CatCont::catIndexString(k) + "]";
					string catMuName = "catMu_part" + catIstr;

					DependentParameter catMuDep(catMuName, catMuGroup, catMuSource);
					gibbs.addParameter(catMuDep);
				}

			}
			else {
				// catMu individual
				MH_Parameter catMu;
				catMu.name = "catMu_part" + catIstr;
				catMu.group = catMuGroup;

				//catMu.deviateFunction = bind(normalDeviate, _1, catMuCandidateSd);
				catMu.deviateFunction = bind(CatCont::MemFirst::catMuDeviateFunction, _1, mhTuningSd.at("catMu_part"), shouldModuloCatMu);
				catMu.llFunction = bind(&CatFirstModel::_LLS_PR_catMu, this, _1, _2, i, k);

				gibbs.addParameter(catMu, catMuStartGrid[k]);
			}

			// ----------------------------------------------
			// catActive
			double catActiveStartValue = (CatCont::uniformDeviate(0, 1) > 0.5) ? 1.0 : 0.0;
			if (this->_catActiveShared) {
				if (i == sharedHolderIndex) {
					MH_Parameter sharedCatActive;
					sharedCatActive.name = "catActive_part" + catIstr;
					sharedCatActive.group = catActiveGroup;

					sharedCatActive.deviateFunction = &CatCont::MemFirst::catActiveDeviateFunction;
					sharedCatActive.llFunction = bind(&CatFirstModel::_LLS_PR_catActive_shared, this, _1, _2, k);

					gibbs.addParameter(sharedCatActive, catActiveStartValue); // Start all categories inactive. TODO
				}
				else {
					string catActiveSource = "catActive_part[" + data.participants.at(sharedHolderIndex).pnum + "," + CatCont::catIndexString(k) + "]";
					string catActiveName = "catActive_part" + catIstr;

					DependentParameter catActiveDep(catActiveName, catActiveGroup, catActiveSource);
					gibbs.addParameter(catActiveDep);
				}
			}
			else {
				// catActive individual
				MH_Parameter catActive;
				catActive.name = "catActive_part" + catIstr;
				catActive.group = catActiveGroup;

				catActive.deviateFunction = &CatCont::MemFirst::catActiveDeviateFunction;
				catActive.llFunction = bind(&CatFirstModel::_LLS_PR_catActive, this, _1, _2, i, k);

				gibbs.addParameter(catActive, catActiveStartValue); //Start all categories inactive. TODO
			}
			
			// ----------------------------------------------
			// catType
			if (this->_catTypeNames.size() == 1) {
				// If only 1 catType, make constant parameter
				ConstantParameter constCatType(catTypeCSIndex, "catType_part" + catIstr, "catType_part");
				gibbs.addParameter(constCatType);

			} 
			else if (this->_catTypeShared) {

				if (i == sharedHolderIndex) {
					MH_Parameter catTypeShared;
					catTypeShared.name = "catType_part" + catIstr;
					catTypeShared.group = "catType_part";

					catTypeShared.deviateFunction = bind(CatFirstModel::_catTypeDeviate, _1, _catTypeNames.size());
					catTypeShared.llFunction = bind(&CatFirstModel::_LLS_PR_catType_shared, this, _1, _2, k);

					gibbs.addParameter(catTypeShared, catTypeCSIndex); // Start catTypes at the cornerstone
				}
				else {
					// If not the holder, catMu and catActive are dependent parameters
					string catTypeSource = "catType_part[" + data.participants.at(sharedHolderIndex).pnum + "," + CatCont::catIndexString(k) + "]";
					string catTypeName = "catType_part" + catIstr;

					DependentParameter catTypeDep(catTypeName, "catType_part", catTypeSource);
					gibbs.addParameter(catTypeDep);
				}

			}
			else {
				// catType individual
				MH_Parameter catType;
				catType.name = "catType_part" + catIstr;
				catType.group = "catType_part";

				catType.deviateFunction = bind(CatFirstModel::_catTypeDeviate, _1, _catTypeNames.size());
				catType.llFunction = bind(&CatFirstModel::_LLS_PR_catType, this, _1, _2, i, k);

				gibbs.addParameter(catType, catTypeCSIndex); // Start catTypes at the cornerstone
			}

		} // k
	} // i

	gibbs.sectionTracker.sectionEnd("Participant Category");
	////////////////////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////////////////////////
	// Condition Effects
	gibbs.sectionTracker.sectionStart("Condition Effects");
	if (runConfig.verbose) {
		CatCont::message("-- Param: Condition Effects");
	}

	vector<string> allParamNames = config.modelConfig.getParamWithAndWithoutConditionEffects();
	vector<string> paramWithCE = config.modelConfig.getParamWithConditionEffects(); // TODO: Instead of this, shouldn't you use the spec DF?
	// Like, no spec results in no CE (constant=0)?

	const map<string, vector<string>>& condSpec = config.cfConfig.condEffects;
	const vector<string>& conditionsForCESpec = condSpec.at("cond"); // Or get from factors?

	for (const string& baseParamName : allParamNames) {

		//CatCont::message("---- Parameter " + baseParamName); // DEBUG

		bool noCE = std::find(paramWithCE.begin(), paramWithCE.end(), baseParamName) == paramWithCE.end();
		bool noClassSpec = condSpec.find(baseParamName) == condSpec.end();

		vector<string> specStr;
		if (noCE || noClassSpec) {
			// DEBUG
			/*
			if (noCE) {
				CatCont::message("---- No condition effects for " + baseParamName);
			}
			if (noClassSpec) {
				CatCont::message("---- No specification found for " + baseParamName);
			}
			*/

			// If not found, constant 0
			// TODO: This should probably be done earlier, like in R
			specStr.resize(config.modelConfig.maxCategories, "constant=0");
		}
		else {
			specStr = condSpec.at(baseParamName);
		}

		//CatCont::message("---- Working out parameter class spec."); // DEBUG

		ParameterClassSpecification pcs;
		pcs.setup(specStr);

		map<string, vector<size_t>> classIndMap = pcs.classIndices();

		for (size_t condInd = 0; condInd < conditionsForCESpec.size(); condInd++) {

			//CatCont::message("---- Condition " + CatCont::toString(condInd)); // DEBUG
			
			string thisParamName = baseParamName + "_cond[" + conditionsForCESpec[condInd] + "]";
			string thisParamGroup = baseParamName + "_cond";

			if (pcs.spec[condInd] == "constant") {

				//CatCont::message("---- Adding constant parameter: " + thisParamName); // DEBUG

				ConstantParameter cp(pcs.numericValue(condInd), thisParamName, thisParamGroup);

				this->gibbs.addParameter(cp); // starting value is ignored

			}
			else if (pcs.spec[condInd] == "class") {
				// Need to find which categories have same class as this.
				// Select the first of them to estimate.

				string thisClassName = pcs.value[condInd];
				const vector<size_t>& thisClassInd = classIndMap.at(thisClassName);
				size_t classEstimatedInd = thisClassInd.front(); // For this class, the condInd that is estimated

				if (classEstimatedInd == condInd) {

					//CatCont::message("---- Finding shared indices"); // DEBUG

					// Get condition indices used by the data for the conditions in this class
					vector<size_t> dataCondInd;
					// For each condition in this class, match it to the condition in data
					for (size_t jClass : thisClassInd) {
						for (size_t jData = 0; jData < this->data.conditionNames.size(); jData++) {
							if (this->data.conditionNames[jData] == conditionsForCESpec[jClass]) {
								dataCondInd.push_back(jData);
								break; // next class index
							}
						}
					}

					//CatCont::message("---- Adding MH parameter: " + thisParamName); // DEBUG

					MH_Parameter mhp;
					mhp.name = thisParamName;
					mhp.group = thisParamGroup;
					mhp.deviateFunction = std::bind(CatCont::normalDeviate, _1, mhTuningSd.at(thisParamGroup));
					mhp.llFunction = std::bind(&CatFirstModel::_LLS_PR_multiCond, this, _1, _2, dataCondInd, baseParamName);

					this->gibbs.addParameter(mhp, 0); // starting values?
				}
				else {

					//CatCont::message("---- Adding dependent parameter: " + thisParamName); // DEBUG

					DependentParameter dp;
					dp.name = thisParamName;
					dp.group = thisParamGroup;
					dp.sourceParameter = baseParamName + "_cond[" + conditionsForCESpec[classEstimatedInd] + "]";
					this->gibbs.addParameter(dp, 0); // starting values?
				}
			}

		}

	}

	gibbs.sectionTracker.sectionEnd("Condition Effects");
	//////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Category Effects (V2, as function of catType.)
	gibbs.sectionTracker.sectionStart("Category Type Effects");
	if (runConfig.verbose) {
		CatCont::message("-- Param: Category Type Effects");
	}

	// TODO: These parameters should have _CTE suffix instead of _cat.

	for (const auto& iter : config.cfConfig.catTypes) {

		if (iter.first == "type") {
			// Skip the "type" column
			continue;
		}

		string baseParamName = iter.first;
		const vector<string>& specStr = iter.second;

		ParameterClassSpecification pcs;
		pcs.setup(specStr);

		map<string, vector<size_t>> classIndMap = pcs.classIndices();

		for (size_t ct = 0; ct < this->_catTypeNames.size(); ct++) {

			// The index for these parameters is the catType name, not the class name.
			string thisParamName = baseParamName + "_cat[" + this->_catTypeNames.at(ct) + "]";
			string thisParamGroup = baseParamName + "_cat";

			if (pcs.spec[ct] == "constant") {

				ConstantParameter cp(pcs.numericValue(ct), thisParamName, thisParamGroup);

				this->gibbs.addParameter(cp); // starting value is ignored

			}
			else if (pcs.spec[ct] == "class") {
				// Find which catTypes have same class as this catType. 
				// Estimate the value of the first of those, the rest being dependent on the first one.

				string thisClassName = pcs.value[ct];
				size_t classEstimatedIndex = classIndMap[thisClassName].front(); // For this class, the catType index that is estimated.

				if (classEstimatedIndex == ct) {
					MH_Parameter mhp;
					mhp.name = thisParamName;
					mhp.group = thisParamGroup;
					mhp.deviateFunction = std::bind(CatCont::normalDeviate, _1, mhTuningSd.at(thisParamGroup));
					mhp.llFunction = std::bind(&CatFirstModel::_LLS_PR_catTypeEffect, this, _1, _2, baseParamName);

					this->gibbs.addParameter(mhp, 0); // starting values?
				}
				else {
					DependentParameter dp;
					dp.name = thisParamName;
					dp.group = thisParamGroup;
					dp.sourceParameter = baseParamName + "_cat[" + this->_catTypeNames.at(classEstimatedIndex) + "]";
					this->gibbs.addParameter(dp, 0); // starting values?
				}
			}
		}
	}
	gibbs.sectionTracker.sectionEnd("Category Type Effects");
	//////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Category Effects (V1, as function of k, not function of catType.)
	/*
	gibbs.sectionTracker.sectionStart("Category Effects");
	if (runConfig.verbose) {
		CatCont::message("-- Param: Category Effects");
	}

	for (const auto& iter : config.cfConfig.catEffects) {

		string baseParamName = iter.first;
		const vector<string>& specStr = iter.second;

		ParameterClassSpecification pcs;
		pcs.setup(specStr);

		map<string, vector<size_t>> classIndMap = pcs.classIndices();

		for (size_t k = 0; k < config.modelConfig.maxCategories; k++) {

			string thisParamName = baseParamName + "_cat[" + CatCont::catIndexString(k) + "]";
			string thisParamGroup = baseParamName + "_cat";

			if (pcs.spec[k] == "constant") {
				
				ConstantParameter cp(pcs.numericValue(k), thisParamName, thisParamGroup);

				this->gibbs.addParameter(cp); // starting value is ignored

			}
			else if (pcs.spec[k] == "class") {
				// Need to find which categories have same class as this.
				// Select the first of them to estimate.

				string thisClassName = pcs.value[k];
				size_t classEstimatedK = classIndMap[thisClassName].front(); // For this class, the k that is estimated

				if (classEstimatedK == k) {
					MH_Parameter mhp;
					mhp.name = thisParamName;
					mhp.group = thisParamGroup;
					mhp.deviateFunction = std::bind(CatCont::normalDeviate, _1, mhTuningSd.at(thisParamGroup));
					mhp.llFunction = std::bind(&CatFirstModel::_LLS_PR_catClass, this, _1, _2, baseParamName);

					this->gibbs.addParameter(mhp, 0); // starting values?
				}
				else {
					DependentParameter dp;
					dp.name = thisParamName;
					dp.group = thisParamGroup;
					dp.sourceParameter = baseParamName + "_cat[" + CatCont::catIndexString(classEstimatedK) + "]";
					this->gibbs.addParameter(dp, 0); // starting values?
				}
			}
		}
	}
	gibbs.sectionTracker.sectionEnd("Category Effects");
	*/
	//////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////
	// Special: Currently just constant values
	gibbs.sectionTracker.sectionStart("Special Parameters");
	if (runConfig.verbose) {
		CatCont::message("-- Param: Special");
	}

	for (size_t i = 0; i < data.participants.size(); i++) {

		string pnum = data.participants.at(i).pnum;

		string name = "catSD_guessShift_part[" + pnum + "]";
		gibbs.addParameter(ConstantParameter(0.0, name, "special"));

		for (size_t k = 0; k < config.modelConfig.maxCategories; k++) {
			string name = "pCatUsedToGuess_part[" + pnum + "," + CatCont::catIndexString(k) + "]";
			gibbs.addParameter(ConstantParameter(1.0, name, "special"));

		}
	}
	gibbs.sectionTracker.sectionEnd("Special Parameters");
	//////////////////////////////////////////////////////////////////////////////

	// TODO: Inappropriate participant parameter likelihoods?
	// Or maybe the new batch likelihood functions could be used instead.

	gibbs.sectionTracker.sectionEnd("Total");

	if (!this->runConfig.profileParameterTypes) {
		gibbs.sectionTracker.clear();
	}
	
	_createParam_modelVariantParameterRemoval();
	
	_doStartingValueOverrides();
	_doConstantValueOverrides();

}


void CatFirstModel::_doConstantValueOverrides(void) {

	const map<string, double> cvo = config.modelConfig.overrides.constantValues;

	if (this->runConfig.verbose) {
		std::stringstream ss;
		ss << "Number of constant parameter value overrides: " << cvo.size();
		CatCont::message(ss.str());
	}

	for (auto it : cvo) {
	//for (const auto& it = cvo.begin(); it != cvo.end(); it++) {
		if (gibbs.hasParameter(it.first)) {

			GibbsParameter* original = gibbs.getParameter<GibbsParameter>(it.first);

			ConstantParameter replacement(it.second, original->name, original->group);

			gibbs.replaceParameter(original->name, replacement);

		} else {
			std::stringstream ss;
			ss << "Note: Constant value provided for \"" << it.first << "\", but there is no parameter by that name.";
			CatCont::message(ss.str(), "_doConstantParameterOverrides");
		}
	}

}

void CatFirstModel::_doStartingValueOverrides(void) {

	const map<string, double>& svo = config.modelConfig.overrides.startingValues;

	if (this->runConfig.verbose) {
		stringstream ss;
		ss << "Number of starting parameter value overrides: " << svo.size();
		CatCont::message(ss.str());
	}

	for (auto it : svo) {
	//for (const auto& it = svo.begin(); it != svo.end(); it++) {
		if (gibbs.hasParameter(it.first)) {

			GibbsParameter* par = gibbs.getParameter<GibbsParameter>(it.first);

			vector<double>& samples = par->getSamples();
			samples.clear();
			samples.push_back(svo.at(par->name));

		}
		else {
			std::stringstream ss;
			ss << "Note: Starting value provided for \"" << it.first << "\", but there is no parameter by that name.";
			CatCont::message(ss.str(), "_doStartingValueOverrides");
		}
	}

}

// Constant value overrides for unused parameters
void CatFirstModel::_createParam_modelVariantParameterRemoval(void) {

	std::map<std::string, double> nameValue;

	// 100 latent is tranformed to 1 manifest, -100 is transformed to 0

	if (config.modelConfig.modelVariant == ModelVariant::BetweenItem) {
		nameValue["pBetween"] = 100; // 1
		nameValue["pContWithin"] = 100; // 1
	}
	else if (config.modelConfig.modelVariant == ModelVariant::WithinItem) {
		nameValue["pBetween"] = -100; // 0
		nameValue["pContBetween"] = 100; // 1
	}
	else if (config.modelConfig.modelVariant == ModelVariant::ZL) {

		//ZL model variant: pBetween = 1, pContBetween = 1, pContWithin = 1, catActive = 0, maxCategories = 0

		nameValue["pBetween"] = 100; //whatever: doesn't have to be 1
		nameValue["pContWithin"] = 100;
		nameValue["pContBetween"] = 100;

		nameValue["pCatGuess"] = -100;

		nameValue["catSelectivity"] = 2; //whatever: doesn't have to be 2
		nameValue["catSD"] = 2;
	}

	for (auto iter : nameValue) {
	//for (auto iter = nameValue.begin(); iter != nameValue.end(); iter++) {
		
		const string& baseParamName = iter.first;
		const double& paramValue = iter.second;

		// Participant parameters
		gibbs.setParameterGroupToConstantValue(baseParamName + "_part", paramValue);

		// Participant population
		string muName = baseParamName + "_part.mu";
		string varName = baseParamName + "_part.var";
		gibbs.replaceParameter(muName, ConstantParameter(0, muName, "populationParam"));
		gibbs.replaceParameter(varName, ConstantParameter(1, varName, "populationParam"));

		// Condition effects
		gibbs.setParameterGroupToConstantValue(baseParamName + "_cond", 0);

		// Category effects
		gibbs.setParameterGroupToConstantValue(baseParamName + "_cat", 0);

	}

}

////////////////////////////////////////////////
// Likelihood functions
// LLS: Log Likelihood Sum (across observations and participants/conditions).
// LLS_PR: LLS with PRior.

double CatFirstModel::_LLS_partByCond(const ParamContainer& param, size_t partInd, size_t condInd) const {

	const ConditionData& condData = this->data.participants.at(partInd).condData.at(condInd);
	if (condData.study.size() == 0) {
		return 0;
	}

	MPMP_complete manifestParam = _paramExtractor.getCompleteMPMP(param, partInd, condInd);

	//MPMP_complete manifestParam;
	//manifestParam.setFromContainer(config.modelConfig, this->_paramExtractor, param, partInd, condInd);

	return _LLS_manifest(condData, manifestParam);
}

double CatFirstModel::_LLS_manifest(const ConditionData& condData, const MPMP_complete& manifestParam) const {
	std::vector<double> likes = _likeCalc.likelihood(condData, manifestParam);

	double lls = 0;
	for (const double& l : likes) {
		lls += std::log(l);
	}
	return lls;
}


double CatFirstModel::_LLS_singleParticipant(const ParamContainer& param, size_t partInd) const {

	const ParticipantData& pData = data.participants.at(partInd);

	double llSum = 0;

	for (size_t j = 0; j < pData.condData.size(); j++) {
		llSum += _LLS_partByCond(param, partInd, j);
	}

	return llSum;
}

double CatFirstModel::_LLS_singleCond(const ParamContainer& param, size_t condInd) const {
	double llSum = 0;

	for (size_t i = 0; i < this->data.participants.size(); i++) {
		llSum += _LLS_partByCond(param, i, condInd);
	}

	return llSum;
}

double CatFirstModel::_LLS_allData(const ParamContainer& param) const {

	double llSum = 0;

	for (size_t i = 0; i < this->data.participants.size(); i++) {
		const ParticipantData& pData = data.participants.at(i);
		for (size_t j = 0; j < pData.condData.size(); j++) {
			llSum += _LLS_partByCond(param, i, j);
		}
	}

	return llSum;
}


////////////////////////////////////////////////////////////////
// LLS_PR functions

double CatFirstModel::_LLS_PR_participantStandard(double thisValue, const ParamContainer& param, size_t partInd, string baseParamName) const {

	double ll = _LLS_singleParticipant(param, partInd);

	double prior = CatCont::normalLL(thisValue, param.at(baseParamName + "_part.mu"), param.at(baseParamName + "_part.var"));

	return ll + prior;

}

double CatFirstModel::_LLS_PR_participantStandard_shared(double thisValue, const ParamContainer& param, double mu0, double var0) const {

	double ll = _LLS_allData(param);

	double prior = CatCont::normalLL(thisValue, mu0, var0);

	return ll + prior;
}

double CatFirstModel::_LLS_PR_multiCond(double thisValue, const ParamContainer& param, vector<size_t> condInds, string baseParamName) const {
	double llSum = 0;
	// Only need to calculate likelihood for conditions that this condition effect parameter is related to.
	for (size_t ci : condInds) {
		llSum += _LLS_singleCond(param, ci);
	}

	//There is only 1 parameter, so only 1 prior.
	double prior = CatCont::cauchyLL(thisValue, priors.at(baseParamName + "_cond.loc"), priors.at(baseParamName + "_cond.scale"));

	return llSum + prior;
}

double CatFirstModel::_LLS_PR_catClass(double thisValue, const ParamContainer& param, string baseParamName) const {

	// No way to split the likelihood
	double llSum = _LLS_allData(param);

	//There is only 1 parameter, so only 1 prior.
	double prior = CatCont::cauchyLL(thisValue, priors.at(baseParamName + "_cat.loc"), priors.at(baseParamName + "_cat.scale"));

	return llSum + prior;
}

double CatFirstModel::_LLS_PR_catTypeEffect(double thisValue, const ParamContainer& param, string baseParamName) const {
	// No way to split the likelihood
	double llSum = _LLS_allData(param);

	//There is only 1 parameter, so only 1 prior.
	double prior = CatCont::cauchyLL(thisValue, priors.at(baseParamName + "_cat.loc"), priors.at(baseParamName + "_cat.scale"));

	return llSum + prior;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Category Parameters (catMu and catActive)

double CatFirstModel::_LLS_PR_catActive(double catActive, const ParamContainer& param, size_t partInd, size_t catIndex) const {

	double ll = _LLS_singleParticipant(param, partInd);

	double prior = _distancePrior_catActive(param, partInd, catIndex);

	// DEBUG
	//string catActiveName = "catActive_part[" + data.participants.at(partInd).pnum + "," + CatCont::catIndexString(catIndex) + "]";
	//CatCont::message(catActiveName + " LL/PR = " + toString(ll) + "/" + toString(prior));

	return ll + prior;
}

double CatFirstModel::_LLS_PR_catMu(double catMu, const ParamContainer& param, size_t partInd, size_t catIndex) const {

	//const ParamContainer::ParamProxy& index = _paramExtractor.part.at(partInd).cat.catActive.at(catIndex);
	//bool catIsActive = param.get(index) == 1.0;

	string catActiveName = "catActive_part[" + data.participants.at(partInd).pnum + "," + CatCont::catIndexString(catIndex) + "]";
	bool catIsActive = param.at(catActiveName) == 1.0;

	// Likelihood and prior only need to be calculated if the category is active
	double llspr = 0;
	if (catIsActive) {
		double ll = _LLS_singleParticipant(param, partInd);

		// If inactive, Although the prior eventually gets calculated as a uniform density, it's faster to just not calculate it.
		double prior = _distancePrior_catMu(param, partInd, catIndex);

		llspr = ll + prior;
	}

	return llspr;
}

// catActive can only be shared if catMu is also shared.
double CatFirstModel::_LLS_PR_catActive_shared(double catActive, const ParamContainer& param, size_t catIndex) const {

	double ll = _LLS_allData(param);

	// If catActive is shared, catMu must be shared as well.
	// pInd=0 to just calculate for shared holder.
	double prior = _distancePrior_catActive(param, 0, catIndex);

	return ll + prior;
}

double CatFirstModel::_LLS_PR_catMu_shared(double catMu, const ParamContainer& param, size_t catIndex) const {

	// likelihood for all participants and conditions
	double ll = this->_LLS_allData(param);

	double prior = 0;
	if (this->_catActiveShared) {
		// If catActive is shared, then catMu and catActive are the same for all participants.
		// Only calculate the prior once, for the holder participant (index 0).
		prior = _distancePrior_catMu(param, 0, catIndex);
	}
	else {
		// If catActive is not shared, the prior is different for each participant.

		// Option 1: Average likelihood. 
		// A category will repel other categories based on the proportion of participants for whom that category is active.
		/*
		double partPriorSum = 0;
		for (size_t i = 0; i < this->data.participants.size(); i++) {
			partPriorSum += _distancePrior_catMu(param, i, catIndex);
		}
		prior = partPriorSum / this->data.participants.size();
		*/

		// Option 2: Do an OR of catActive for each k.
		// If the category is active for any participant, it will repel other categories.
		vector<unsigned int> catActiveUnion(config.modelConfig.maxCategories, 1);
		for (size_t i = 0; i < this->data.participants.size(); i++) {
			CategoryParam cp = _paramExtractor.getParticipantCategory(param, i);
			for (size_t k = 0; k < config.modelConfig.maxCategories; k++) {
				catActiveUnion[k] *= cp.catActive[k];
			}
		}
		CategoryParam cp = _paramExtractor.getParticipantCategory(param, 0);
		prior = _catPriorCalc.distancePrior_catMu(cp.catMu, catActiveUnion, catIndex, true); // log = true

	}

	return ll + prior;
}

double CatFirstModel::_distancePrior_catMu(const ParamContainer& param, size_t partInd, size_t catIndex) const {

	CategoryParam cp = _paramExtractor.getParticipantCategory(param, partInd);
	
	return _catPriorCalc.distancePrior_catMu(cp.catMu, cp.catActive, catIndex, true); // log = true
}

double CatFirstModel::_distancePrior_catActive(const ParamContainer& param, size_t partInd, size_t catIndex) const {

	CategoryParam cp = _paramExtractor.getParticipantCategory(param, partInd);

	// DEBUG
	//CatCont::message("catMu/catActive size = " + CatCont::toString(cp.catMu.size()) + "/" + CatCont::toString(cp.catActive.size()));
	//CatCont::message("Active catMu size = " + CatCont::toString(cp.activeCatMu.size()));

	return _catPriorCalc.distancePrior_catActive(cp.catMu, cp.catActive, catIndex, true); // log = true
}



double CatFirstModel::_LLS_PR_catType(double catType, const ParamContainer& param, size_t partInd, size_t catIndex) const {

	string catActiveName = "catActive_part[" + data.participants.at(partInd).pnum + "," + CatCont::catIndexString(catIndex) + "]";
	bool catIsActive = param.at(catActiveName) == 1.0;

	// Likelihood only needs to be calculated if the category is active
	double ll = 0;
	if (catIsActive) {
		ll = _LLS_singleParticipant(param, partInd);
	}

	// Look up prior based on catType, which is treated as an index
	double prior = this->priors.at("catType_part[" + _catTypeNames.at((size_t)catType) + "]");
	prior = log(prior);
	// TODO: These probabilities could be copied into a vector and looked up based on index for a small speed increase.
	// Then would not need to save _catTypeNames.

	return ll + prior;
}

double CatFirstModel::_LLS_PR_catType_shared(double catType, const ParamContainer& param, size_t catIndex) const {

	double ll = 0;

	if (this->_catActiveShared) {
		// Likelihood only needs to be calculated if the category is active.
		string catActiveName = "catActive_part[" + data.participants.at(0).pnum + "," + CatCont::catIndexString(catIndex) + "]";
		bool catIsActive = param.at(catActiveName) == 1.0;
		if (catIsActive) {
			ll = _LLS_allData(param);
		}
	}
	else {
		// Only need to calculate likelihood for participants for whom this category is active.
		for (size_t i = 0; i < this->data.participants.size(); i++) {

			string catActiveName = "catActive_part[" + data.participants.at(i).pnum + "," + CatCont::catIndexString(catIndex) + "]";
			bool catIsActive = param.at(catActiveName) == 1.0;

			if (catIsActive) {
				ll += _LLS_singleParticipant(param, i);
			}
		}
	}

	// Look up prior based on catType, which is treated as an index
	double prior = this->priors.at("catType_part[" + _catTypeNames.at((size_t)catType) + "]");
	prior = log(prior);

	return ll + prior;
}

// TODO: Move to somewhere (CCM_Util?)
size_t sampleIntExclusive(size_t exclude, size_t nValues) {
	// If nValues == 4, sample a value from [0,2.999]
	double unifMax = nValues - (1 + 1e-12);
	double samp = R::runif(0, unifMax);

	size_t downgrade = floor(samp); // Convert sample to integer by taking its floor

	// If at the excluded value or higher
	if (downgrade >= exclude) {
		downgrade++; // Step past the excluded value
	}
	return downgrade;
}

double CatFirstModel::_catTypeDeviate(double currentType, size_t typeCount) {
	// currentType must be <= typeCount - 1
	return sampleIntExclusive((size_t)currentType, typeCount);
}

///////////////////////////////////////////////////////////////
// Conjugate sampling functions for hyperparameters on latent participant parameters.

double CatFirstModel::_conjugateSample_mu(const ParamContainer& param, string baseParamName, double mu0, double var0) const {
	vector<double> y = this->gibbs.getCurrentGroupValues(baseParamName + "_part");
	double currentVar = param.at(baseParamName + "_part.var");
	return normal_muPostSample(y, currentVar, mu0, var0);
}

double CatFirstModel::_conjugateSample_var(const ParamContainer& param, string baseParamName, double a0, double b0) const {
	vector<double> y = this->gibbs.getCurrentGroupValues(baseParamName + "_part");
	double currentMu = param.at(baseParamName + "_part.mu");
	return normal_varPostSample(y, currentMu, a0, b0);
}


/////////////////////////////////////////////////////////////
// Priors and MH Tuning

map<string, double> CatFirstModel::getDefaultPriors(const vector<string>& catTypeNames) {

	map<string, double> defPriors;

	vector<string> probParams = { "pMem", "pBetween", "pContBetween", "pContWithin", "pCatGuess" };

	for (const string& pp : probParams) {

		// Hierarchical priors on participant parameters
		defPriors[pp + "_part.mu.mu"] = 0;
		defPriors[pp + "_part.mu.var"] = pow(1.2, 2); //sd = 1.2, var = 1.44
		defPriors[pp + "_part.var.a"] = 1.5;
		defPriors[pp + "_part.var.b"] = 1.5;

		// Priors on condition effects
		defPriors[pp + "_cond.loc"] = 0;
		defPriors[pp + "_cond.scale"] = 0.3;

		// Priors on category effects
		defPriors[pp + "_cat.loc"] = 0;
		defPriors[pp + "_cat.scale"] = 0.3;

	}

	// TODO: Apply linear range to SD params and catMu

	vector<string> sdParams = { "contSD", "catSD", "catSelectivity" };

	for (const string& sp : sdParams) {

		// Hierarchical priors on participant parameters
		defPriors[sp + "_part.mu.mu"] = 35;
		defPriors[sp + "_part.mu.var"] = pow(15, 2); //sd = 15
		defPriors[sp + "_part.var.a"] = 0.75;
		defPriors[sp + "_part.var.b"] = 0.75;

		// Priors on condition effects
		defPriors[sp + "_cond.loc"] = 0;
		defPriors[sp + "_cond.scale"] = 3;

		// Priors on category effects
		defPriors[sp + "_cat.loc"] = 0;
		defPriors[sp + "_cat.scale"] = 3;

	}

	defPriors["catMuPriorSD"] = 12;
	//defPriors["catActivePriorProb"] = 0.5;

	// These don't favor the cornerstone, they are just placeholders
	for (const string& ctn : catTypeNames) {
		defPriors["catType_part[" + ctn + "]"] = 1.0 / catTypeNames.size();
	}

	return defPriors;
}

map<string, double> CatFirstModel::getDefaultMHTuning(void) {

	// TODO: Apply linear range to SD params and catMu for both participant and condition params.

	// These are standard deviations of normal candidate distributions.
	map<string, double> defTuningSD;

	///////////////////////////////////////////
	// Participant parameters. All are latent.
	defTuningSD["pMem_part"] = 0.3;
	defTuningSD["pBetween_part"] = 0.9;
	defTuningSD["pContBetween_part"] = 0.4;
	defTuningSD["pContWithin_part"] = 0.7;
	defTuningSD["pCatGuess_part"] = 0.7;

	// SD parameters and catMu are in degrees
	defTuningSD["contSD_part"] = 2;
	defTuningSD["catSD_part"] = 1;
	defTuningSD["catSelectivity_part"] = 2;

	defTuningSD["catMu_part"] = 5;

	///////////////////////////////////////////
	// Condition effects. All are latent.
	defTuningSD["pMem_cond"] = 0.2;
	defTuningSD["pBetween_cond"] = 0.4;
	defTuningSD["pContBetween_cond"] = 0.2;
	defTuningSD["pContWithin_cond"] = 0.2;
	defTuningSD["pCatGuess_cond"] = 0.3;

	defTuningSD["contSD_cond"] = 1;
	defTuningSD["catSD_cond"] = 0.5;
	defTuningSD["catSelectivity_cond"] = 0.8;

	// Category effects. Latent
	defTuningSD["pMem_cat"] = 0.2;
	defTuningSD["pBetween_cat"] = 0.4;
	defTuningSD["pContBetween_cat"] = 0.2;
	defTuningSD["pContWithin_cat"] = 0.2;
	//defTuningSD["pCatGuess_cond"] = 0.3;

	defTuningSD["contSD_cat"] = 1;
	defTuningSD["catSD_cat"] = 0.5;
	//defTuningSD["catSelectivity_cond"] = 0.8;

	return defTuningSD;
}

void CatFirstModel::_setPriors(void) {
	
	this->priors = CatFirstModel::getDefaultPriors(config.cfConfig.catTypes.at("type"));

	// Do prior overrides
	if (this->runConfig.verbose) {
		stringstream ss;
		ss << "Number of prior overrides: " << config.modelConfig.overrides.priors.size();
		CatCont::message(ss.str());
	}

	for (const auto& it : config.modelConfig.overrides.priors) {
	//for (auto it = config.overrides.priors.begin(); it != config.overrides.priors.end(); it++) {

		string name = it.first;

		if (this->priors.find(name) != this->priors.end()) {
			this->priors[name] = it.second;
		} else {
			CatCont::message("Warning: Invalid prior override key: \"" + name + "\" ignored.", "_setPriors");
		}

	}

	//After the overrides are in, calculate these things:
	_catPriorCalc.setup(config.modelConfig, this->priors["catMuPriorSD"]);
}

void CatFirstModel::_setMHTuning(void) {

	this->mhTuningSd = CatFirstModel::getDefaultMHTuning();

	// Do MH overrides
	if (this->runConfig.verbose) {
		stringstream ss;
		ss << "Number of MH overrides: " << config.modelConfig.overrides.mhTunings.size();
		CatCont::message(ss.str());
	}

	for (auto it : config.modelConfig.overrides.mhTunings) {
	//for (auto it = config.overrides.mhTunings.begin(); it != config.overrides.mhTunings.end(); it++) {

		const string& name = it.first;

		if (mhTuningSd.find(name) != mhTuningSd.end()) {
			mhTuningSd[name] = it.second;
		} else {
			CatCont::message("Warning: Invalid MH tuning override key: \"" + name + "\" ignored.", "_setMHTuning");
		}
	}

}



} // namespace CatFirst
} // namespace CatCont
