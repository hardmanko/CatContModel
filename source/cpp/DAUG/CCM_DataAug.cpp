#include "CCM_DataAug.h"

namespace CatCont {

	void daugModel::_createParameters(void) {

		using namespace std::placeholders;

		gibbs.clear();

		///////////////////////////////////////
		// Observation Labels

		// Every observation needs a label: pnum,cond,obs
		vector<string> obsIndex(3);

		for (size_t i = 0; i < data.participants.size(); i++) {
			const ParticipantData& thisPart = data.participants.at(i);

			obsIndex[0] = thisPart.pnum;

			for (size_t j = 0; j < thisPart.condData.size(); j++) {

				obsIndex[1] = thisPart.condData.at(j).condition;

				for (size_t obs = 0; obs < thisPart.condData.at(j).response.size(); obs++) {

					obsIndex[2] = _convertToString(obs + 1); // 1 indexed

					double startLabel = 0; // TODO

					{
						MH_Parameter obsLabel;
						obsLabel.name = joinIndexedParam("label", catIndex);
						obsLabel.group = "label";

						obsLabel.llFunction = std::bind(&daugModel::_obsLL_MH, this, _1, _2, i, j, obs);
						obsLabel.deviateFunction = std::bind(&daugModel::_labelDeviate, this, _1);

						gibbs.addParameter(obsLabel, startLabel);
					}


				}
			}
		}

		///////////////////////////
		// Label Proportions


		

	}

	double _singleObservationLL(string label) const {

		string labelBase = extractBaseParameterName(label);

		CombinedParameters param;
		if (labelBase == "mem_cont") {
			param.pMem = 1;
			param.pContBetween = 0;
			param.pContWithin = 0; // ???
			param.pCatGuess = 0; // doesn't matter really
		} else if (labelBase == "guess_unif") {
			param.pMem = 0;
			param.pCatGuess = 0;
		} else if (labelBase == "mem_cat") {
			param.pMem = 1;
			param.pContBetween = 1;
			param.pContWithin = 1; // ????
		}
		
		// The issue with reusing the KVR17 likelihood is that all categories are needed to calculate the weights,
		// but only the labeled category distribution should be used to calulate the daug likelihood.

		//vector<double> likelihood = Circular::betweenAndWithinLikelihood(param, conditionData, modelConfig);

	}

	// Instead of index, can use just response as arg?
	double daugModel::_obsLL_MH(double newLabel, const ParameterList& param, size_t partIndex, size_t condIndex, size_t obsIndex) const {
		double response = data.participants.at(partIndex).condData.at(condIndex).response.at(obsIndex);

		// TODO

		return 1;
	}

	double daugModel::_labelDeviate(double currentLabel) const {

		double max = dcfg.indexedLabels.size() - 0.00001;

		double newLabel = currentLabel;
		while (newLabel == currentLabel) {
			newLabel = std::floor(R::runif(0, max));
		}

		return newLabel;
	}

	double daugModel::_labelNameToIndex(string name) const {

	}

	string daugModel::_labelIndexToName(double index) const {

	}



} // namespace CatCont