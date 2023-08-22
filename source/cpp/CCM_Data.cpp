#include "CCM_Data.h"

namespace CatCont {


vector<ParticipantData> copyParticipantData(vector<string> pnumsCol, vector<string> condsCol, vector<double> studyCol, vector<double> responseCol, DataType dataType, bool verbose) {

	vector<string> uniquePnums = uniqueElements(pnumsCol);
	vector<string> uniqueConditions = uniqueElements(condsCol);

	vector<size_t> trialsPerCondition(uniqueConditions.size(), 0);

	vector<ParticipantData> data;

	// V1
	for (size_t i = 0; i < uniquePnums.size(); i++) {

		ParticipantData thisPart;
		thisPart.pnum = uniquePnums[i];

		// Find rows of data for this participant
		vector<size_t> pnumRows;
		for (size_t row = 0; row < pnumsCol.size(); row++) {
			if (pnumsCol[row] == uniquePnums[i]) {
				pnumRows.push_back(row);
			}
		}

		for (size_t j = 0; j < uniqueConditions.size(); j++) {

			//Copy from the secondary data frame to primitive data types in a ConditionData struct
			ConditionData condData;

			condData.condition = uniqueConditions[j];

			for (size_t obs = 0; obs < pnumRows.size(); obs++) {

				size_t thisRow = pnumRows[obs];

				//If the condition on the current row is equal to the jth condition, store the data
				if (condsCol[thisRow] == uniqueConditions[j]) {
					double study = studyCol[thisRow];
					double response = responseCol[thisRow];

					if (dataType == DataType::Circular) {
						study = Circular::degreesToRadians(study); //TODO: This should probably be done by the model or in R
						response = Circular::degreesToRadians(response);
					}
					condData.study.push_back(study);
					condData.response.push_back(response);
					condData.originalRow.push_back(thisRow);
				}
			}

			thisPart.condData.push_back(condData);

			trialsPerCondition[j] += condData.study.size();

			/*
			if (verbose) {
				std::stringstream ss;
				ss << "Participant " << uniquePnums[i] << ", condition " << uniqueConditions[j] << ": " <<
					_convertToString(condData.study.size()) << " trials found.";
				CatCont::message(ss.str());
			}
			*/
		}

		data.push_back(thisPart);

	}

	if (verbose) {
		std::stringstream ss;
		ss << "Participants found: " << uniquePnums.size() << endl;
		ss << "Total trials per condition: " << endl;
		for (unsigned int j = 0; j < uniqueConditions.size(); j++) {
			//string cond = uniqueConditions[j];
			ss << uniqueConditions[j] << ": " << trialsPerCondition[j] << endl;
		}
		CatCont::message(ss.str());
	}

	return data;

}

#if COMPILING_WITH_RCPP
bool DataCollection::setup(Rcpp::DataFrame df, DataType dataType, bool verbose) {

	// Copy the columns from R data frame
	Rcpp::CharacterVector pnumsColRaw = df["pnum"];
	vector<string> pnumsCol(pnumsColRaw.begin(), pnumsColRaw.end());

	Rcpp::CharacterVector condsColRaw = df["cond"];
	vector<string> condsCol(condsColRaw.begin(), condsColRaw.end());

	Rcpp::NumericVector studyColRaw = df["study"];
	vector<double> studyCol(studyColRaw.begin(), studyColRaw.end());

	Rcpp::NumericVector respColRaw = df["response"];
	vector<double> respCol(respColRaw.begin(), respColRaw.end());

	const vector<string> uniquePnums = uniqueElements(pnumsCol);
	const vector<string> uniqueConditions = uniqueElements(condsCol);

	vector<size_t> trialsPerCondition(uniqueConditions.size(), 0);

	// Reset this
	*this = DataCollection();

	this->conditionNames = uniqueConditions;

	// TODO: No need to alias here
	//vector<ParticipantData>& partData = this->participants;

	// V2: Loop once through obs

	// Map from string to index in this->participants
	map<string, size_t> pnumIndexMap;
	map<string, size_t> condIndexMap;

	// Initialize participants
	this->participants.resize(uniquePnums.size());
	for (size_t i = 0; i < uniquePnums.size(); i++) {
		pnumIndexMap[uniquePnums[i]] = i;

		this->participants[i].pnum = uniquePnums[i];
		this->participants[i].condData.resize(uniqueConditions.size());
		for (size_t j = 0; j < uniqueConditions.size(); j++) {
			this->participants[i].condData[j].condition = uniqueConditions[j];
		}
	}

	for (size_t j = 0; j < uniqueConditions.size(); j++) {
		condIndexMap[uniqueConditions[j]] = j;
	}

	// Loop once
	for (size_t obs = 0; obs < pnumsCol.size(); obs++) {
		size_t i = pnumIndexMap[pnumsCol[obs]];
		size_t j = condIndexMap[condsCol[obs]];

		ConditionData& condData = this->participants[i].condData[j];

		condData.study.push_back(studyCol[obs]);
		condData.response.push_back(respCol[obs]);
		condData.originalRow.push_back(obs);

		trialsPerCondition[j]++;
	}

	// Once data are in, convert units and set the ranges
	this->degreesToRadians(dataType);

	this->setDataRanges(this->participants);


	if (verbose) {
		std::stringstream ss;
		ss << "Participants found: " << uniquePnums.size() << endl;
		ss << "Total trials per condition: " << endl;
		for (size_t j = 0; j < uniqueConditions.size(); j++) {
			ss << uniqueConditions[j] << ": " << trialsPerCondition[j] << endl;
		}
		CatCont::message(ss.str());
	}

	return true;
}

void DataCollection::degreesToRadians(DataType dataType) {

	if (dataType != DataType::Circular) {
		return;
	}

	for (size_t i = 0; i < this->participants.size(); i++) {
		for (size_t j = 0; j < this->conditionNames.size(); j++) {
			this->participants[i].condData[j].study = Circular::degreesToRadians(this->participants[i].condData[j].study);
			this->participants[i].condData[j].response = Circular::degreesToRadians(this->participants[i].condData[j].response);
		}
	}
}

// This is the same as setup(), but uses copyParticipantData
bool DataCollection::setFromDataFrame(Rcpp::DataFrame df, DataType dataType, bool verbose) {

	// Copy the columns from R data frame
	Rcpp::CharacterVector pnumsColRaw = df["pnum"];
	vector<string> pnumsCol(pnumsColRaw.begin(), pnumsColRaw.end());

	Rcpp::CharacterVector condsColRaw = df["cond"];
	vector<string> condsCol(condsColRaw.begin(), condsColRaw.end());

	Rcpp::NumericVector studyColRaw = df["study"];
	vector<double> studyCol(studyColRaw.begin(), studyColRaw.end());

	Rcpp::NumericVector respColRaw = df["response"];
	vector<double> respCol(respColRaw.begin(), respColRaw.end());

	// Reset this
	*this = DataCollection();

	this->participants = copyParticipantData(pnumsCol, condsCol, studyCol, respCol, dataType, verbose);

	// Set condition names
	this->conditionNames = uniqueElements(condsCol);

	this->setDataRanges(this->participants);

	return true;
}
#endif

void DataCollection::setDataRanges(const vector<ParticipantData>& partData) {
	for (size_t i = 0; i < this->participants.size(); i++) {
		for (size_t j = 0; j < this->conditionNames.size(); j++) {

			const ConditionData& cd = this->participants[i].condData[j];

			for (size_t obs = 0; obs < cd.study.size(); obs++) {

				this->studyRange.lower = min(this->studyRange.lower, cd.study[obs]);
				this->studyRange.upper = max(this->studyRange.upper, cd.study[obs]);

				this->responseRange.lower = min(this->responseRange.lower, cd.response[obs]);
				this->responseRange.upper = max(this->responseRange.upper, cd.response[obs]);
			} // obs
		} // j
	} // i
}


} // namespace CatCont