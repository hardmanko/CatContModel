#pragma once

#include <set>

#include "CCM_Main.h"
#include "CCM_Data.h"
#include "CCM_Util.h"
#include "CCM_Circular.h"

#include "R_Interface.h"

#include "CF_ModelConfig.h"
#include "CF_Parameters.h"
#include "CF_Likelihood.h"

#include "MF_ModelUtil.h"

Rcpp::List CCM_CF_RPE(Rcpp::DataFrame dataDF, Rcpp::List modelConfig, Rcpp::List runConfig);

namespace CatCont {
namespace CatFirst {


class CatFirstModel {
public:

	bool setup(Rcpp::DataFrame dataDF, Rcpp::List modelConfig, Rcpp::List runConfig);

	//bool RPE(Rcpp::DataFrame dataDF, Rcpp::List modelConfig, Rcpp::List runConfig, Rcpp::List equalityConstraints);


	//ModelConfiguration modelConfig;
	//CatFirstConfig CFConfig;
	//VariablePrecisionConfig VPConfig;
	CF_CompositeConfig config;

	RunConfig runConfig;

	DataCollection data;


	GibbsSampler gibbs;
	void createGibbsParameters(void);

	void _createParam_modelVariantParameterRemoval(void);
	void _doConstantValueOverrides(void);
	void _doStartingValueOverrides(void);


	// Maybe just implement these in R?
	static map<string, double> getDefaultPriors(const vector<string>& catTypeNames);
	static map<string, double> getDefaultMHTuning(void);

	map<string, double> priors; //fixed priors
	map<string, double> mhTuningSd;

	void _setPriors(void);
	void _setMHTuning(void);


private:

	ParamExtractor _paramExtractor;
	LikelihoodCalculator _likeCalc;

	// Log Likelihood Sum functions
	double _LLS_partByCond(const ParamContainer& param, size_t partInd, size_t condInd) const;
	double _LLS_manifest(const ConditionData& condData, const MPMP_complete& manifestParam) const;
	
	double _LLS_singleParticipant(const ParamContainer& param, size_t partInd) const;
	double _LLS_singleCond(const ParamContainer& param, size_t condInd) const;
	double _LLS_allData(const ParamContainer& param) const;

	// LLS with PRiors
	double _LLS_PR_participantStandard(double thisValue, const ParamContainer& param, size_t partInd, string baseParamName) const;
	double _LLS_PR_participantStandard_shared(double thisValue, const ParamContainer& param, double mu0, double var0) const;

	double _LLS_PR_multiCond(double thisValue, const ParamContainer& param, vector<size_t> condInd, string baseParamName) const;

	double _LLS_PR_catClass(double thisValue, const ParamContainer& param, string baseParamName) const;
	double _LLS_PR_catTypeEffect(double thisValue, const ParamContainer& param, string baseParamName) const;

	// Category parameters
	MemFirst::CategoryParamPriorCalculator _catPriorCalc;
	double _distancePrior_catMu(const ParamContainer& param, size_t partInd, size_t catIndex) const;
	double _distancePrior_catActive(const ParamContainer& param, size_t partInd, size_t catIndex) const;
	//bool _catMuShared = false;
	//bool _catActiveShared = false;

	double _LLS_PR_catMu(double catMu, const ParamContainer& param, size_t partInd, size_t catIndex) const;
	double _LLS_PR_catMu_shared(double catMu, const ParamContainer& param, size_t catIndex) const;
	double _LLS_PR_catActive(double catActive, const ParamContainer& param, size_t partInd, size_t catIndex) const;
	double _LLS_PR_catActive_shared(double catActive, const ParamContainer& param, size_t catIndex) const;

	bool _catActiveShared = false;
	bool _catTypeShared = false;
	vector<string> _catTypeNames; // TODO: Where does this thing go?
	double _LLS_PR_catType(double catType, const ParamContainer& param, size_t partInd, size_t catIndex) const;
	double _LLS_PR_catType_shared(double catType, const ParamContainer& param, size_t catIndex) const;
	static double _catTypeDeviate(double currentType, size_t typeCount);


	// Conjugate updating functions for param_part.mu and .var
	double _conjugateSample_mu(const ParamContainer& param, string baseParamName, double mu0, double var0) const;
	double _conjugateSample_var(const ParamContainer& param, string baseParamName, double a0, double b0) const;


};

} // namespace CatFirst
} // namespace CatCont
