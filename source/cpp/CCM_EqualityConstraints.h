#pragma once

#include "CCM_Util.h"

namespace CatCont {

//Condition parameters are the natural use for this.
//Some participant parameters could be constrained with this: catSD = contSD, for example.
struct EqualityConstraints {

    typedef map<string, string> IndividualMappings;
    typedef map<string, map<unsigned int, vector<unsigned int>>> GroupMappings;

    static string FreeParameter;

    //make private?
    IndividualMappings individualMappings;
    GroupMappings groupMappings;

    bool setup(const IndividualMappings& mappings, const vector<string>& conditionNames, const vector<string>& allParameterNames);


    const vector<unsigned int>& getEqualConditionIndices(string param, unsigned int cond) const;
    
	bool isFreeParameter(string parameter) const;
    string getSourceParameter(string parameter) const;

    IndividualMappings incorporateAdditionalParameters(IndividualMappings mappings, const vector<string>& allNames) const;

    bool simplifyIndividualMappings(IndividualMappings* mapping) const;

    GroupMappings calculateEqualConditionIndices(const IndividualMappings& mapping, const vector<string>& conditionNames) const;
};

struct ParameterClassSpecification {

	// TODO: Cornerstone is strange and should maybe be done in R.
	//static const CornerstoneClassName = "CS";

	// trackingIndex is like cond name or cat index. Parallels the specifications.
	bool setup(vector<string> specifications) {

		/*
		if (specifications.size() != trackingIndex_.size()) {
			CatCont::stop("Length of specifications and trackingIndex do not match.", "ParameterClassSpecification::setup");
		}

		this->trackingIndex = trackingIndex_;
		*/

		this->spec.resize(specifications.size());
		this->value.resize(specifications.size());

		for (size_t i = 0; i < specifications.size(); i++) {
			vector<string> split = CatCont::splitString(specifications[i], "=");

			// Other checks are all done in R
			if (split.size() != 2) {
				CatCont::stop("Invalid specification: " + specifications[i], "ParameterClassSpecification::setup");
			}
			this->spec[i] = split[0];
			this->value[i] = split[1];
		}

		return true;
	}

	//vector<string> trackingIndex;
	vector<string> spec;
	vector<string> value;

	double numericValue(size_t index) const {
		return CatCont::fromString<double>(value.at(index));
	}

	vector<string> uniqueClasses(void) const {
		vector<string> classValues;
		for (size_t i = 0; i < this->spec.size(); i++) {
			if (this->spec[i] == "class") {
				classValues.push_back(this->value[i]);
			}
		}
		return CatCont::uniqueElements(classValues);
	}

	//vector<size_t> indicesOfSameClass(string className) const; // maybe?

	// map<class_name, class_indices>
	map<string, vector<size_t>> classIndices(void) const {
		map<string, vector<size_t>> rval;
		for (size_t i = 0; i < this->spec.size(); i++) {
			if (this->spec[i] == "class") {
				rval[this->value[i]].push_back(i);
			}
		}
		return rval;
	}

};

} // namespace CatCont