#include "CCM_EqualityConstraints.h"

namespace CatCont {
    
string EqualityConstraints::FreeParameter = "FREE_PARAMETER";

bool EqualityConstraints::setup(const IndividualMappings& mappings, const vector<string>& conditionNames, const vector<string>& allParameterNames) {
    individualMappings = mappings;

    if (allParameterNames.size() > 0) {
        vector<string> keysToErase;

        for (auto it = individualMappings.begin(); it != individualMappings.end(); it++) {

            bool firstNotFound = find(allParameterNames.begin(), allParameterNames.end(), it->first) == allParameterNames.end();
            bool secondNotFound = find(allParameterNames.begin(), allParameterNames.end(), it->second) == allParameterNames.end();

            if (it->first == FreeParameter) {
                firstNotFound = false;
            }
            if (it->second == FreeParameter) {
                secondNotFound = false;
            }

            if (firstNotFound) {
                CatCont::message("Parameter " + it->first + " does not exist and will be ignored.", "EqualityConstraints");
            }

            if (secondNotFound) {
                CatCont::message("Parameter " + it->second + " does not exist and will be ignored.", "EqualityConstraints");
            }

            if (firstNotFound || secondNotFound) {
                keysToErase.push_back(it->first);
            }
        }

        for (string s : keysToErase) {
            individualMappings.erase(s);
        }
    }

    bool simplifySuccess = simplifyIndividualMappings(&individualMappings);
    if (!simplifySuccess) {
        return false;
    }

    groupMappings = calculateEqualConditionIndices(individualMappings, conditionNames);
    return true;
}



//This should include the given condition
const vector<unsigned int>& EqualityConstraints::getEqualConditionIndices(string param, unsigned int cond) const {
    return groupMappings.at(param).at(cond);
}

bool EqualityConstraints::isFreeParameter(string parameter) const {
	return getSourceParameter(parameter) == FreeParameter;
}

//Basically, read out of individualMappings
string EqualityConstraints::getSourceParameter(string parameter) const {
    if (individualMappings.find(parameter) == individualMappings.end()) {
        return "PARAMETER_NOT_FOUND";
    }

    return individualMappings.at(parameter);
}

//Account for parameters not yet mentioned(?)
EqualityConstraints::IndividualMappings EqualityConstraints::incorporateAdditionalParameters(IndividualMappings mappings, const vector<string>& allNames) const {
    for (string s : allNames) {
        if (mappings.find(s) == mappings.end()) {
            mappings[s] = FreeParameter;
        }
    }
    return mappings;
}

//Reduce to simplest form
bool EqualityConstraints::simplifyIndividualMappings(IndividualMappings* mapping) const {

    //get parameter names
    vector<string> names;
    for (auto it = mapping->begin(); it != mapping->end(); it++) {
        names.push_back(it->first);
    }

    bool mappingChanged = false;
    for (const string& currentTarget : names) {

        string currentSource = mapping->at(currentTarget);

        if (currentSource == FreeParameter) {
            continue;
        }

        set<string> cycleHistory;
        cycleHistory.insert(currentSource);

        string mostSimpleSource = currentSource;
        while (true) {
            string nextPossibleSource = mapping->at(mostSimpleSource);

            if (cycleHistory.find(nextPossibleSource) != cycleHistory.end()) {
                CatCont::message("Error: Equality constraint cycle detected.", "EqualityConstraints");
                return false;
            } else {
                cycleHistory.insert(nextPossibleSource);
            }

            if (nextPossibleSource != FreeParameter) {
                mostSimpleSource = nextPossibleSource;
            } else {
                break;
            }
        }
        if (mostSimpleSource != currentSource) {
            mapping->at(currentTarget) = mostSimpleSource;
            mappingChanged = true;
        }

    }

    if (mappingChanged) {
        return simplifyIndividualMappings(mapping);
    }

    return true;
}



//conditionNames comes from the DataCollection struct.
EqualityConstraints::GroupMappings EqualityConstraints::calculateEqualConditionIndices(const IndividualMappings& mapping, const vector<string>& conditionNames) const {

    GroupMappings rval;

    map<string, unsigned int> cnameToIndex;
    for (unsigned int i = 0; i < conditionNames.size(); i++) {
        cnameToIndex[conditionNames[i]] = i;
    }

    vector<string> names;
    for (auto it = mapping.begin(); it != mapping.end(); it++) {
		// Only take parameters with _cond.
		// TODO: Shouldn't it be that free parameters are rejected?
        if (it->first.find("_cond[") != string::npos) {
            names.push_back(it->first);
        }
    }

    for (const string& possibleSource : names) {

        string baseParameterName = extractBaseParameterName(possibleSource);

        string sourceCondName = extractIndex(possibleSource);

        vector<string> cNames;
        vector<unsigned int> cInds;

        //If this is a free parameter, add it to its own group
        if (mapping.at(possibleSource) == FreeParameter) {
            string sourceCondName = extractIndex(possibleSource);

            cNames.push_back(sourceCondName);
            cInds.push_back(cnameToIndex.at(sourceCondName));
        }

        for (const string& possibleTarget : names) {

            if (mapping.at(possibleTarget) == possibleSource) {
                string targetCondName = extractIndex(possibleTarget);

                cNames.push_back(targetCondName);
                cInds.push_back(cnameToIndex.at(targetCondName));
            }

        }

        unsigned int sourceConditionIndex = cnameToIndex.at(sourceCondName);
        rval[baseParameterName][sourceConditionIndex] = cInds;
    }

    return rval;
}

} // namespace CatCont