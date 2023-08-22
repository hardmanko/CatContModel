#include "MF_Parameters.h"

using namespace std;

namespace CatCont {
namespace MemFirst {


ParticipantParameters getParticipantParameters(const ParameterList& param, string pnum, unsigned int maxCategories) {
	string istr = "_part[" + pnum + "]";

	ParticipantParameters part;
	part.pMem = param.at("pMem" + istr);
	part.pBetween = param.at("pBetween" + istr);
	part.pContBetween = param.at("pContBetween" + istr);
	part.pContWithin = param.at("pContWithin" + istr);
	part.pCatGuess = param.at("pCatGuess" + istr);

	part.contSD = param.at("contSD" + istr);

	string catIStrBase = "_part[" + pnum + ",";

	for (unsigned int k = 0; k < maxCategories; k++) {
		string catIstr = catIStrBase + CatCont::catIndexString(k) + "]";

		double active = param.at("catActive" + catIstr);

		if (active == 1.0) {
			part.cat.mu.push_back(param.at("catMu" + catIstr));
		}
		//if not active, don't add it to the list
	}

	part.cat.selectivity = param.at("catSelectivity" + istr);
	part.cat.SD = param.at("catSD" + istr);

	return part;
}

ConditionParameters getConditionParameters(const ParameterList& param, string condName) {

	ConditionParameters cond;

	string istr = "_cond[" + condName + "]";

	cond.pMem = param.at("pMem" + istr);
	cond.pBetween = param.at("pBetween" + istr);
	cond.pContBetween = param.at("pContBetween" + istr);
	cond.pContWithin = param.at("pContWithin" + istr);
	cond.pCatGuess = param.at("pCatGuess" + istr);

	cond.contSD = param.at("contSD" + istr);
	cond.cat.selectivity = param.at("catSelectivity" + istr);
	cond.cat.SD = param.at("catSD" + istr);

	return cond;
}


ParticipantParameters getParticipantParameters(const ParamContainer& param, string pnum, unsigned int maxCategories) {
	string istr = "_part[" + pnum + "]";

	ParticipantParameters part;
	part.pMem = param.at("pMem" + istr);
	part.pBetween = param.at("pBetween" + istr);
	part.pContBetween = param.at("pContBetween" + istr);
	part.pContWithin = param.at("pContWithin" + istr);
	part.pCatGuess = param.at("pCatGuess" + istr);

	part.contSD = param.at("contSD" + istr);

	string catIStrBase = "_part[" + pnum + ",";

	for (unsigned int k = 0; k < maxCategories; k++) {
		string catIstr = catIStrBase + CatCont::catIndexString(k) + "]";

		double active = param.at("catActive" + catIstr);

		if (active == 1.0) {
			part.cat.mu.push_back(param.at("catMu" + catIstr));
		}
		//if not active, don't add it to the list
	}

	part.cat.selectivity = param.at("catSelectivity" + istr);
	part.cat.SD = param.at("catSD" + istr);

	return part;
}

ConditionParameters getConditionParameters(const ParamContainer& param, string condName) {
	ConditionParameters cond;

	string istr = "_cond[" + condName + "]";

	cond.pMem = param.at("pMem" + istr);
	cond.pBetween = param.at("pBetween" + istr);
	cond.pContBetween = param.at("pContBetween" + istr);
	cond.pContWithin = param.at("pContWithin" + istr);
	cond.pCatGuess = param.at("pCatGuess" + istr);

	cond.contSD = param.at("contSD" + istr);
	cond.cat.selectivity = param.at("catSelectivity" + istr);
	cond.cat.SD = param.at("catSD" + istr);

	return cond;
}

// Combine participant and condition effects.
// Transform parameters from latent to manifest spaces: Probs to [0,1], SDs from SD degrees to precision radians, and catMu from degrees to radians.
CombinedParameters combineParameters(const ParticipantParameters& part, const ConditionParameters& cond, const SDRanges& sdRanges, DataType dataType) {
	CombinedParameters rval;

	rval.pMem = CatCont::paramTransform_probability(part.pMem + cond.pMem);
	rval.pBetween = CatCont::paramTransform_probability(part.pBetween + cond.pBetween);
	rval.pContBetween = CatCont::paramTransform_probability(part.pContBetween + cond.pContBetween);
	rval.pContWithin = CatCont::paramTransform_probability(part.pContWithin + cond.pContWithin);
	rval.pCatGuess = CatCont::paramTransform_probability(part.pCatGuess + cond.pCatGuess);

	rval.contSD = CatCont::paramTransform_sd(part.contSD + cond.contSD, sdRanges, dataType);
	rval.cat.selectivity = CatCont::paramTransform_sd(part.cat.selectivity + cond.cat.selectivity, sdRanges, dataType);
	rval.cat.SD = CatCont::paramTransform_sd(part.cat.SD + cond.cat.SD, sdRanges, dataType);

	// catMu are estimated in radians. No condition parameters or transformations are applied. 
	// The only possible transformation would be fmod(mu, 2*PI), but it is unnecessary, so skipped to increase estimation speed.
	//rval.cat.mu = part.cat.mu;

	// Estimate catMu in degrees rather than radians. (Or linear units)
	// This saves on the need to constantly convert back and forth between R and C++.
	rval.cat.mu = CatCont::Circular::degreesToRadians(part.cat.mu);
	//rval.cat.mu = CatCont::Circular::clampAngle360(rval.cat.mu); // Not needed. Candidate catMu are clamped to 360

	return rval;
}

} // namespace MemFirst
} // namespace CatCont