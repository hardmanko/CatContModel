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
    
    string getSourceParameter(string parameter) const;

    IndividualMappings incorporateAdditionalParameters(IndividualMappings mappings, const vector<string>& allNames) const;

    bool simplifyIndividualMappings(IndividualMappings* mapping) const;

    GroupMappings calculateEqualConditionIndices(const IndividualMappings& mapping, const vector<string>& conditionNames) const;
};

} // namespace CatCont