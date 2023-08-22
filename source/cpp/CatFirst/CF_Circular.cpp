#include "CF_Model.h"


namespace CatCont {
namespace CatFirst {

vector<CompleteParam> calculateParamPerK(bool activeOnly, const ModelConfiguration& modelConfig, const ParticipantParam& part, const StandardParam& cond, const CompleteCatClassParam& catClass) {
	
	vector<size_t> calcK;
	if (activeOnly) {
		calcK = part.cat.getActiveCatK();
	} else {
		for (size_t i = 0; i < modelConfig.maxCategories; i++) {
			calcK.push_back(i);
		}
	}

	vector<CompleteParam> paramPerK(modelConfig.maxCategories);
	for (size_t i = 0; i < calcK.size(); i++) {
		CompleteParam cp;
		cp.cat = part.cat;

		cp.special = part.special; // TODO: Is this at all right? Where does special come from?
		// TODO: Special param modifications to standard?

		cp.standard = part.standard + cond + catClass.param[k];

		cp.standard.transformToManifest(modelConfig.ranges, modelConfig.dataType);

		paramPerK.push_back(cp);
	}

}

vector<double> CatFirstModel::likelihood_raw_singleIter(const ModelConfiguration& modelConfig, const ConditionData& data, const ParameterList& parList, string pnum, string cond)
{

	ParticipantParam part;
	part.standard.setFromList(parList, "_part[" + pnum + "]");
	part.cat.setFromList(parList, pnum, modelConfig.maxCategories);
	//part.special = TODO

	StandardParam cond;
	cond.setFromList(parList, "_cond[" + cond + "]");

	CompleteCatClassParam ccp;
	ccp.setFromList(parList, ccm);


}

// Should only be calculated once based on parameters. Does not depend on data.
// rval is length ak
vector<double> calculateCatGuessWeights(const MPMP_guessing& guessPar) {

	size_t nCatActive = guessPar.cat.nCatActive();

	vector<double> catGuessWeights(nCatActive);

	double catGuessWeight_sum = 0;

	for (size_t ak = 0; ak < nCatActive; ak++) {

		size_t k = guessPar.cat.activeCatInd[ak];

		catGuessWeights[ak] = guessPar.pCatUsedToGuess[k] / nCatActive;
		catGuessWeight_sum += catGuessWeights[ak];
	}
	// Normalize catGuessWeights
	for (size_t ak = 0; ak < nCatActive; ak++) {
		catGuessWeights[ak] /= catGuessWeight_sum;
	}

	return catGuessWeights;
}

// catGuess_catWeight is length ak
double likelihood_guessing(double response, const MPMP_guessing& guessPar, const vector<double>& catGuessWeights) {
	// Calculate guessing component of the model, which does not depend on memory
	double guessDens_weightedSum = 0;
	for (size_t ak = 0; ak < nCatActive; ak++) {

		size_t k = catParam.activeCatInd[ak];

		//const CompleteParam& cp = paramPerK[k];

		//double catSD_guess = cp.standard.catSD + cp.special.catSD_guessShift;
		double catGuessDens = vmLut.dVonMises(response, catParam.catMu[k], guessPar.catSD[k]);

		double guessDens = guessPar.pCatGuess * catGuessDens + (1 - guessPar.pCatGuess) * unifDens;

		double weightedGuessDens = catGuessWeights[ak] * guessDens;

		guessDens_weightedSum += weightedGuessDens;
	}
	return guessDens_weightedSum;
}

double likelihood_memory_BI_singleK(double study, double response, const MPMP_memory& par, const vector<double>& studyWeights) {



}

vector<double> CatFirstModel::likelihood_BI(const ModelConfiguration& modelConfig, const ConditionData& data,
	const ParticipantParam& part, const StandardParam& cond, const CompleteCatClassParam& catClass)
{

	const double zeroWeightThreshold = 1e-10; // move to somewhere

	const double unifDens = 1 / (2 * PI); // move to constexpr function or define

	// Alias categorization param
	const CategorizationParam& catParam = part.cat; // TODO: From where
	size_t nCatActive = catParam.nCatActive();

	// Prepare category weights vector (which is recalculated for each study value)

	vector<double> studyWeights(nCatActive);
	//vector<double> responseWeights(nCatActive);

	vector<double> likelihoods(data.study.size());
	//double weightedLikeSum = 0;

	// TODO: When guessing, what parameter values are used? Which K are they drawn from?
	// pCatGuess does not depend on k.
	MPMP_guessing guessingParam;
	guessingParam.setCategories(catParam);
	//guessingParam.pCatGuess = ? ? ? ? ;

	vector<MPMP_memory> memParam_perAK(nCatActive);
	for (size_t ak = 0; ak < nCatActive; ak++) {
		MPMP_memory mem;
		
		StandardParam combinedStandard = part.standard + cond + catClass.param[k];
		combinedStandard.transformToManifest(modelConfig.ranges, modelConfig.dataType);

		mem.pMem = combinedStandard.pMem;
		mem.pBetween = combinedStandard.pBetween;
		mem.pContBetween = combinedStandard.pContBetween;
		mem.pContWithin = combinedStandard.pContWithin;
		mem.contSD = combinedStandard.contSD;
		mem.catSD = combinedStandard.catSD;

		// combined for condition effects
		mem.catSelectivity = combinedStandard.catSelectivity;

		mem.cat = part.cat; // TODO: not part.cat, but maybe some other way?

		// TODO: More?

		memParam_perAK[ak] = mem;

		// Set guessing parameters
		guessingParam.catSD[ak] = mem.catSD + part.special.catSD_guessShift; // TODO: Where do special come from?
		//guessingParam.pCatUsedToGuess[ak] = ? ? ? ? ; 
	}

	// TODO: This only needs to be done for active categories, which is a small savings.
	vector<CompleteParam> paramPerK;
	for (size_t k = 0; k < modelConfig.maxCategories; k++) {
		CompleteParam cp;
		cp.cat = part.cat; // copy

		cp.special = part.special; // TODO: Is this at all right? Where does special come from?
		// TODO: Special param modifications to standard?

		cp.standard = part.standard + cond + catClass.param[k];

		cp.standard.transformToManifest(modelConfig.ranges, modelConfig.dataType);

		paramPerK.push_back(cp);
	}

	//vector<double> activeCats = 

	// Calculate the probability (weight) that each category will be used when making a catGuess.
	// This depends on using the bizarro pCatUsedToGuess parameter.
	vector<double> catGuess_catWeight(nCatActive);
	double catGuess_catWeight_sum = 0;
	for (size_t ak = 0; ak < nCatActive; ak++) {

		size_t k = catParam.activeCatInd[ak];

		catGuess_catWeight[ak] = cp.special.pCatUsedToGuess[k] / nCatActive;
		catGuess_catWeight_sum += catGuess_catWeight[ak];
	}
	// Normalize catGuess_catWeight
	for (size_t ak = 0; ak < nCatActive; ak++) {
		catGuess_catWeight[ak] /= catGuess_catWeight_sum;
	}


	// Calculate likelihood for each observation.
	for (size_t obs = 0; obs < data.study.size(); obs++) {

		// Calculate guessing component of the model, which does not depend on memory
		double guessDens_weightedSum = 0;
		for (size_t ak = 0; ak < nCatActive; ak++) {

			size_t k = catParam.activeCatInd[ak];

			const CompleteParam& cp = paramPerK[k];

			double catSD_guess = cp.standard.catSD + cp.special.catSD_guessShift;
			double catGuessDens = vmLut.dVonMises(data.response[obs], catParam.catMu[k], catSD_guess);

			double guessDens = cp.standard.pCatGuess * catGuessDens + (1 - cp.standard.pCatGuess) * unifDens;

			double weightedGuessDens = catGuess_catWeight[ak] * guessDens;

			guessDens_weightedSum += weightedGuessDens;
		}



		// Calculate categorization weights for this observation
		// TODO: Weights function needs correct args.
		CatCont::Circular::categoryWeights(data.study[obs], par, studyWeights.data());
		//CatCont::Circular::categoryWeights(data.response[obs], par, responseWeights.data());

		// Calculate density within each categorization (unless weight is 0, in which case skip)
		vector<double> catLikelihood(nCatActive);
		double obsDens_weightedSum = 0;
		if (nCatActive == 0) {
			// If no cats are active, calculate like ZL


		}
		else {
			for (size_t ak = 0; ak < nCatActive; ak++) {

				// If the study weight is 0, then the following tree has probability 0
				if (studyWeights[ak] < zeroWeightThreshold) {
					catLikelihood[ak] = 0;
					continue;
				}

				// k is the original category index
				size_t k = catParam.activeCatInd[ak];

				const CompleteParam& cp = paramPerK[k];

				// Density for categorial memory response. Centered on category.
				double catMemDens = vmLut.dVonMises(data.response[obs], catParam.activeCatMu[ak], cp.standard.catSD);

				// Density for continuous memory response. Centered on study vlue.
				double contMemDens = vmLut.dVonMises(data.response[obs], data.study[obs], cp.standard.contSD);

				double memoryDensity = (cp.standard.pContBetween * contMemDens) + ((1 - cp.standard.pContBetween) * catMemDens);

				// Add guessing (which is the same for all cats, but pMem may be different)
				double unweightedDens = cp.standard.pMem * memoryDensity + (1 - cp.standard.pMem) * guessDens_weightedSum;

				catLikelihood[ak] = studyWeights[ak] * unweightedDens;

				obsDens_weightedSum += catLikelihood[ak];
			} // ak
		}

		//The PI / 180 scales the likelihood so that it would be correct if degrees had been used instead of radians
		likelihoods[obs] = obsDens_weightedSum * (PI / 180.0);

	}

	return likelihoods;

} // likelihood_BI()

static vector<double> CatFirstModel::likelihood_WI(const ModelConfiguration& modelConfig, const ConditionData& data, const MPMP_complete& par) {

	const double unifDens = 1 / (2 * PI); // move to constexpr function or define

	size_t nCatActive = par.guessing.cat->nCatActive(); // TODO: Maybe have one more global copy of cat?

	// Precalculate catGuessWeights because they do not depend on data.
	vector<double> catGuessWeights = calculateCatGuessWeights(par.guessing);

	// Precalculate withinKappas because they do not depend on the data.
	vector<double> withinKappa_ak(nCatActive);
	for (size_t i = 0; i < nCatActive; i++) {
		if (config.withinItem_contSD_is_withinSD) {
			// In this special case, withinKappa is based on only contSD, so contSD is withinSD
			withinKappa_ak[i] = par.memory[i].contSD; // contSD has already been converted into precision. (Bad variable naming)

		}
		else {
			// Normallly, contSD and catSD are combined based on pContWithin.
			withinKappa_ak[i] = Circular::combineKappas(par.memory[i].pContWithin, par.memory[i].contSD, par.memory[i].catSD);
		}
	}

	// When calculating the within component predicted response center, weight the 
	// categorical and continuous components of the response by pContWithin
	vector<double> cmLoc(2);
	vector<double> cmWeights(2);
	cmWeights[0] = par.pContWithin;
	cmWeights[1] = 1 - par.pContWithin;


	// Calculate likelihood for each observation.
	vector<double> likelihoods(data.study.size());
	for (size_t obs = 0; obs < data.study.size(); obs++) {

		double L_guess = likelihood_guessing(data.response[obs], par.guessing, catGuessWeights);


	}

	return likelihoods;
}

} // namespace CatFirst
} // namespace CatCont