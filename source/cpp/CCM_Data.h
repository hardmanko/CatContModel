// Things related to reading data

#pragma once

#include "CCM_Main.h"

#include "CCM_Util.h"
#include "CCM_ModelConfig.h"

#include "CCM_Circular.h"

namespace CatCont {

	struct ConditionData {
		string condition; // redundant with DataCollection::conditionNames

		vector<double> study;
		vector<double> response;
		vector<size_t> originalRow; // TODO: One way to track data source.
	};

	struct ParticipantData {
		string pnum;
		vector<ConditionData> condData;
	};

	struct DataCollection {
		DataCollection(void) {
			studyRange.lower = numeric_limits<double>::max();
			studyRange.upper = numeric_limits<double>::min();

			responseRange.lower = numeric_limits<double>::max();
			responseRange.upper = numeric_limits<double>::min();
		}

		struct {
			double lower;
			double upper;
		} studyRange;

		struct {
			double lower;
			double upper;
		} responseRange;

		void setDataRanges(const vector<ParticipantData>& partData);
		void degreesToRadians(DataType dataType);

		vector<ParticipantData> participants;
		vector<string> conditionNames;

#if COMPILING_WITH_RCPP
		bool setup(Rcpp::DataFrame df, DataType dataType, bool verbose);
		bool setFromDataFrame(Rcpp::DataFrame df, DataType dataType, bool verbose);

#endif
#if COMPILING_WITH_CX
		bool setFromDelimFile(string filename, DataType dataType, bool verbose);
#endif
	};

	vector<ParticipantData> copyParticipantData(vector<string> pnumsCol, vector<string> condsCol,
		vector<double> studyCol, vector<double> responseCol, DataType dataType, bool verbose);

	/*
	#if COMPILING_WITH_CX
		vector<ParticipantData> getParticipantData(string filename, DataType dataType, bool verbose);
	#endif
	#if COMPILING_WITH_RCPP
		vector<ParticipantData> getParticipantData(Rcpp::DataFrame df, DataType dataType, bool verbose);
	#endif
	*/

} // namespace CatCont