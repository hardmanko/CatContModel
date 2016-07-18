#pragma once

#include "GibbsSampler.h"
#include "UtilityFunctions.h"

#include "CCM_Util.h"
#include "CCM_Linear.h"
#include "CCM_Circular.h"

#define VON_MISES_MAX_SD 10000
#define VON_MISES_STEP_SIZE 0.01

namespace CatCont {

	class Bayesian {
	public:

		struct SDRanges {
			double minSd;
			double maxSd;

			double minPrecision;
			double maxPrecision;
		};

		struct Configuration {
			string outputDir; //Not used with Rcpp

			unsigned int iterations;

			DataType dataType;
			ModelVariant modelVariant;

			unsigned int maxCategories;
			string cornerstoneConditionName;
			unsigned int cornerstoneConditionIndex;

			SDRanges ranges;

			bool calculateParticipantLikelihoods;


			vector<string> paramWithConditionEffects;

			vector<string> getParamWithAndWithoutConditionEffects(void) const;
			vector<string> getParamWithConditionEffects(void) const;
			vector<string> getParamWithoutConditionEffects(void) const;
			vector<string> getParamWithHierachicalPriors(void) const;


			Linear::LinearConfiguration linearConfiguration;
		};

		Bayesian();

		void setup(string configFilename);

		void setData(vector<ParticipantData> data);
		void setPriors(void);
		void setMhTuning(void);
		void createParameters(void);

		void run(void);
		void outputResults(void);

		GibbsSampler gibbs;

		Data data;

		map<string, double> mhTuningSd;
		map<string, double> priors; //fixed priors

#ifdef COMPILING_WITH_CX
		map<string, string> configFile;
#endif
#ifdef COMPILING_WITH_RCPP
		struct {
			Rcpp::List mhTuningOverrides;
			Rcpp::List priorOverrides;
			Rcpp::List startingValueOverrides;
			Rcpp::List constantValueOverrides;
		} rcppConfig;
#endif

		Configuration config;

		bool iterationCallback(GibbsSampler* gs, unsigned int iteration);

		void setParameterStartingValues(map<string, double> vals);

		static ParticipantParameters getParticipantParameters(const ParameterList& param, string pnum, unsigned int maxCategories);
		static ConditionParameters getConditionParameters(const ParameterList& param, string condName);
		static CombinedParameters combineParameters(const ParticipantParameters& part, const ConditionParameters& cond, const SDRanges& ranges, DataType dataType);

	private:

		void _doMhOverrides(void);
		void _doPriorOverrides(void);
		void _doStartingValueOverrides(void);
		void _doConstantParameterOverrides(void);


		double singleParticipant_ll(const ParameterList& param, unsigned int pIndex) const;
		double singleCond_ll(const ParameterList& param, unsigned int condIndex) const;

		double normalPriorParameter_ll(double thisParam, const ParameterList& param, unsigned int pIndex, string paramName) const;


		double postMuSample(const ParameterList& param, string paramSetName, double mu0, double var0) const;
		double postVarSample(const ParameterList& param, string paramSetName, double a0, double b0) const;

		double genericConditionParameter_ll(double thisParam, const ParameterList& param, unsigned int condIndex, string baseName) const;

		static double _sdParameterTransformation(double sd, const SDRanges& ranges, DataType dataType);
		static double _catActiveDeviateFunction(double active);
		
		//Convenience versions of the static functions with similar names
		ParticipantParameters _getParticipantParameters(const ParameterList& param, unsigned int pIndex) const;
		ConditionParameters _getConditionParameters(const ParameterList& param, unsigned int condIndex) const;
		CombinedParameters _combineParameters(const ParticipantParameters& part, const ConditionParameters& cond) const;


		double _llFunction(const CombinedParameters& par, const ConditionData& data) const;

		double _singleParticipantLLSamplingFunction(const ParameterList& param, unsigned int i) const;


		struct {
			double sd; //standard deviation in degrees/units
			double kappa; //precision in radians
			double maxLikelihood; //non-log

			unsigned int catMuScaleFactorApproximationSteps;
		} _catMuPriorData;

		double _scaledCatMuDensity(const vector<double>& mus, const vector<unsigned int>& catActives, unsigned int k, unsigned int steps) const;
		double _catMuPenaltyDensity(const vector<double>& mus, const vector<unsigned int>& catActive, unsigned int k) const;
		double _estimateCatMuScaleFactor(vector<double> mus, const vector<unsigned int>& catActives, unsigned int k, unsigned int steps) const;
		double _catMuPenalityPrior(const ParameterList& param, unsigned int pIndex, unsigned int catIndex) const;
		double _catActivePenaltyPrior(const ParameterList& param, unsigned int pIndex, unsigned int catIndex) const;
		
		double catActive_ll(double catActive, const ParameterList& param, unsigned int pIndex, unsigned int catIndex) const;
		double catMu_ll(double catMu, const ParameterList& param, unsigned int pIndex, unsigned int catIndex) const;

		ParameterList _getDecorrelatingValues(const ParameterList& param, string paramSetName, double deviate) const;
		ParameterList decorrelatingCurrent(const ParameterList& param, string paramSetName) const;
		ParameterList decorrelatingCandidate(const ParameterList& param, string paramSetName, double mhSd) const;
		double decorrelating_ll(const ParameterList& theseParam, const ParameterList& allParam, string paramSetName) const;

	};

	map<string, map<string, double>> calculateWAIC(const vector<ParticipantData>& allData, const Bayesian::Configuration& modelConfig, const vector< ParameterList >& posteriorIterations);
}
