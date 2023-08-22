#pragma once

#include "CCM_Main.h"
#include "CCM_ModelConfig.h"
#include "CCM_Data.h"

#include "GibbsSampler.h"
#include "UtilityFunctions.h"

#include "CCM_Util.h"
#include "CCM_Linear.h"
#include "CCM_Circular.h"
#include "CCM_EqualityConstraints.h"
#include "CCM_DistributionLUTs.h"

#include "MF_Parameters.h"
#include "MF_ModelUtil.h"
#include "MF_Likelihood.h"

namespace CatCont {
	namespace MemFirst {

map<string, map<string, double>> calculateWAIC(
	//const vector<ParticipantData>& allData, 
	const DataCollection& allData,
	const ModelConfiguration& modelConfig, 
	const vector< ParamMap >& posteriorIterations);


class MemFirstModelRunner {
public:


	//void setData(vector<ParticipantData> data);

	void createParameters(void);

	DataCollection data;

	RunConfig runConfig;
	ModelConfiguration config;

	map<string, double> mhTuningSd;
	map<string, double> priors; //fixed priors

	GibbsSampler gibbs;

#ifdef COMPILING_WITH_CX
	map<string, string> configFile;
#endif


	//static ParticipantParameters getParticipantParameters(const ParamContainer& param, string pnum, unsigned int maxCategories);
	//static ConditionParameters getConditionParameters(const ParamContainer& param, string condName);
	//static CombinedParameters combineParameters(const ParticipantParameters& part, const ConditionParameters& cond, const SDRanges& sdRanges, DataType dataType);

	//bool iterationCallback(GibbsSampler* gs, unsigned int iteration);

	static map<string, double> getDefaultPriors(void);
	static map<string, double> getDefaultMHTuning(void);

private:

	// Some cached info
	bool _catActiveShared = false;

	// Setup
	void _setPriors(void);
	void _setMhTuning(void);

	void _doMhOverrides(void);
	void _doPriorOverrides(void);
	void _doStartingValueOverrides(void);
	void _doConstantParameterOverrides(void);

	// Likelihood functions and wrappers
	double _llFunction(const CombinedParameters& par, const ConditionData& data) const;

	double _allData_ll(const ParamContainer& param) const;

	double singleParticipant_ll(const ParamContainer& param, unsigned int pIndex) const;
	double singleCond_ll(const ParamContainer& param, unsigned int condIndex) const;

	double normalPriorParameter_ll(double thisParam, const ParamContainer& param, unsigned int pIndex, string paramName) const;
	//double genericConditionParameter_ll(double thisParam, const ParamContainer& param, unsigned int condIndex, string baseName) const;
	double multiConditionParameter_ll(double thisParam, const ParamContainer& param, vector<unsigned int> condIndices, string baseName) const;

	//double _singleParticipantLLSamplingFunction(const ParamContainer& param, unsigned int pIndex) const;

	double postMuSample(const ParamContainer& param, string baseParamName, double mu0, double var0) const;
	double postVarSample(const ParamContainer& param, string baseParamName, double a0, double b0) const;


	//Convenience versions of the static functions with similar names
	ParticipantParameters _getParticipantParameters(const ParamContainer& param, unsigned int pIndex) const;
	ConditionParameters _getConditionParameters(const ParamContainer& param, unsigned int condIndex) const;
	CombinedParameters _combineParameters(const ParticipantParameters& part, const ConditionParameters& cond) const;

	

	// Category parameter functions and data
	/*
	struct {
		double sd; //standard deviation in degrees/units
		double kappa; //precision in radians
		double maxLikelihood; //non-log
	} _catMuPriorData;

	double _scaledCatMuDensity(const vector<double>& mus, const vector<unsigned int>& catActives, unsigned int catIndex, unsigned int steps) const;
	double _catMuPenaltyDensity(const vector<double>& mus, const vector<unsigned int>& catActive, unsigned int catIndex) const;
	double _estimateCatMuScaleFactor(vector<double> mus, const vector<unsigned int>& catActives, unsigned int catIndex, unsigned int steps) const;
	double _catMuPenalityPrior(const ParamContainer& param, unsigned int pIndex, unsigned int catIndex) const;
	double _catActivePenaltyPrior(const ParamContainer& param, unsigned int pIndex, unsigned int catIndex) const;
	*/

	CategoryParamPriorCalculator _catPriorCalc;

	double _distancePrior_catMu(const ParamContainer& param, unsigned int pIndex, unsigned int catIndex) const;
	double _distancePrior_catActive(const ParamContainer& param, unsigned int pIndex, unsigned int catIndex) const;

	double _catActive_ll(double catActive, const ParamContainer& param, unsigned int pIndex, unsigned int catIndex) const;
	double _catMu_ll(double catMu, const ParamContainer& param, unsigned int pIndex, unsigned int catIndex) const;

	double _catMu_shared_ll(double catMu, const ParamContainer& param, unsigned int catIndex) const;
	double _catActive_shared_ll(double catActive, const ParamContainer& param, unsigned int catIndex) const;

	// Decorellating steps to reduce cross-correlations between parameters
	ParamMap _getDecorrelatingValues(const ParamContainer& param, string baseParamName, double deviate) const;
	ParamMap decorrelatingCurrent(const ParamContainer& param, string baseParamName) const;
	ParamMap decorrelatingCandidate(const ParamContainer& param, string baseParamName, double mhSd) const;
	double decorrelating_ll(const ParamMap& theseParam, const ParamContainer& allParam, string baseParamName) const;

};

} // namespace MemFirst
} // namespace CatCont
