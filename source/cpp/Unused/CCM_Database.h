#pragma once

#include <functional>

#include "CCM_Util.h"

namespace CatCont {

	struct StudyResponseData {
		vector<double> study;
		vector<double> response;
	};

	struct PartCondCell {
		StudyResponseData cellData;

		CombinedParameters combinedParam;

		vector< vector< double> > cachedWeights;
	};

	struct PartNode {
		std::string pnumStr;

		ParticipantParameters partParam;

		//vector<PartCondCell> conds;
	};

	struct CondNode {
		std::string condName;

		ConditionParameters condParam;
	};



	enum class ParamInd : int {
		pMem,
		pBetween,
		pContBetween,
		pContWithin,
		pCatGuess,

		contSD,

		catSD,
		catSel,

		catMu,
		catActive,

		psPlatW,
		psSplineW
	};





	/* This class could save processing time by reducing the amount of string-based parameter lookups.

	The gibbs sampler needs to have functionality added to it for this class to work. The gibbs sampler
	would need to have a way to use functions from getPartParamGetter(), for example, to get and set parameter
	values.

	*/
	class CCM_Database {
	public:

		struct Config {

			std::function<double(double)> probParamTrans;
			std::function<double(double)> probParamTransInverse;

			std::function<double(double)> sdParamTrans;
			std::function<double(double)> sdParamTransInverse;

		};

		bool setup(const Config& cfg);
		void addParticipant();

		using PartInd = int; // or size_t or whatever
		using CondInd = int;
		using CatInd = int;

		const CombinedParameters& getCombinedParam(PartInd part, CondInd cond);

		std::function<double(void)> getPartParamGetter(PartInd part, ParamInd param);
		std::function<void(double)> getPartParamSetter(PartInd part, ParamInd param);

		std::function<double(void)> getCatParamGetter(PartInd part, ParamInd param, CatInd k);
		std::function<void(double)> getCatParamSetter(PartInd part, ParamInd param, CatInd k);

		std::function<double(void)> getCondParamGetter(CondInd cond, ParamInd param);
		std::function<void(double)> getCondParamSetter(CondInd cond, ParamInd param);

		std::string getParamStr(ParamInd param, PartInd part);
		std::string getCatParamStr(ParamInd param, PartInd part, CatInd k);

		void setPartParam(PartInd part, ParamInd param, double val);

		void setCatParam(PartInd part, ParamInd param, CatInd k, double val);

	private:

		Config _cfg;

		std::vector<PartNode> _parts;
		std::vector<CondNode> _conds;
		std::vector< std::vector<PartCondCell> > _cells;

		std::map<std::string, PartInd> _partStrToInd;
		std::map<std::string, CondInd> _condStrToInd;

	};

} // namespace CatCont