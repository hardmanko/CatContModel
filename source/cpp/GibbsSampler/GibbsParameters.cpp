#include "GibbsParameters.h"

#include "GibbsSampler.h"

std::uniform_real_distribution<double> GibbsSampler::canonical = std::uniform_real_distribution<double>(0, 1);

////////////////////
// GibbsParameter //
////////////////////

void GibbsParameter::_updateCurrentValue(double s) {
	_gibbs->setCurrentParameterValue(this->name, s);
}

//////////////////
// MH_Parameter //
//////////////////
double MH_Parameter::_getNextSample(void) {

	ParameterList param = _gibbs->getCurrentParameterValues();

	double currentValue = this->value();
	double candidate = this->deviateFunction(currentValue);

	//if the candidate is outside of the range of possible values, reject it
	if (candidate < range.lower || candidate > range.upper) {
		acceptanceTracker.parameterRejected();
		return currentValue;
	}

	param[this->name] = currentValue; // what? why would this->value() not equal the stored value?
	_gibbs->updateDependentParameters(&param); //Does this need to be here?
	double currentLL = this->llFunction(currentValue, param);

	param[this->name] = candidate;
	_gibbs->updateDependentParameters(&param);
	double candidateLL = this->llFunction(candidate, param);


	double likelihoodRatio = exp(candidateLL - currentLL);

	double newValue;
	double canonicalSample = GibbsSampler::canonical(_gibbs->getGenerator());

	//if the likelihood ratio is greater than 1, the sample is always less than it, so accept
	if (canonicalSample < likelihoodRatio) {
		newValue = candidate;
		acceptanceTracker.parameterAccepted();
	} else {
		newValue = currentValue;
		acceptanceTracker.parameterRejected();
	}

	return newValue;
}

////////////////////////
// DependentParameter //
////////////////////////

double DependentParameter::evaluate(const ParameterList& param) const {
	if (sourceParameter == name) {
		return _samples.back();
	}

	return param.at(sourceParameter);

	return this->_gibbs->getParameter<GibbsParameter>(sourceParameter)->value();
}

double DependentParameter::value(void) const {

	//If this is its own source, return the last sample
	if (sourceParameter == name) {
		return _samples.back();
	}

	return this->_gibbs->getParameter<GibbsParameter>(sourceParameter)->value();
}

double DependentParameter::getSample(unsigned int iteration) const {
	return this->_gibbs->getParameter<GibbsParameter>(sourceParameter)->getSample(iteration);
}

std::vector<double>& DependentParameter::getSamples(void) {
	if (sourceParameter == name) {
		return _samples;
	}

	return this->_gibbs->getParameter<GibbsParameter>(sourceParameter)->getSamples();
}

double DependentParameter::_getNextSample(void) {
	return this->value();
}

/////////////////////////
// CalculatedParameter //
/////////////////////////

double CalculatedParameter::_getNextSample(void) {
	return this->value();
}

double CalculatedParameter::value(void) const {
	return evaluate(_gibbs->getCurrentParameterValues());
}

double CalculatedParameter::evaluate(const ParameterList& param) const {
	if (samplingFunction != nullptr) {
		return samplingFunction(param);
	}
	throw(std::runtime_error("CalculatedParameter not configured with a sampling function!"));
	return 0;
}

double CalculatedParameter::getSample(unsigned int iteration) const {
	return samplingFunction(_gibbs->getIterationParameterValues(iteration));
}

std::vector<double>& CalculatedParameter::getSamples(void) {

	//for (unsigned int i = 0; i < _samples.size(); i++) {
	//
	//}

	return _samples;
}


///////////////////////
// DecorrelatingStep //
///////////////////////
double DecorrelatingStep::_getNextSample(void) {

	ParameterList paramCopy = _gibbs->getCurrentParameterValues();

	std::map<std::string, double> currentValues = currentValuesFunction(paramCopy);
	double llCurrent = llFunction(currentValues, paramCopy);

	std::map<std::string, double> candidateValues = candidateValuesFunction(paramCopy);
	for (const auto& p : candidateValues) {
		paramCopy[p.first] = p.second;
	}
	_gibbs->updateDependentParameters(&paramCopy);
	double llCandidate = llFunction(candidateValues, paramCopy);


	double likelihoodRatio = exp(llCandidate - llCurrent);

	//double newValue;
	double canonicalSample = GibbsSampler::canonical(_gibbs->getGenerator());

	//if the likelihood ratio is greater than 1, the sample is always less than it, so accept
	if (canonicalSample < likelihoodRatio) {
		acceptanceTracker.parameterAccepted();

		for (const auto& p : candidateValues) {

			UNSAFE_storeNewSample(_gibbs->getParameter<GibbsParameter>(p.first), p.second);

			//if accepting, change the current values of the parameters stored by the sampler
			//if not accepting, nothing needs to be done to the current values as they have not been changed
			_gibbs->setCurrentParameterValue(p.first, p.second);
		}
	} else {
		acceptanceTracker.parameterRejected();

		for (const auto& p : currentValues) {
			UNSAFE_storeNewSample(_gibbs->getParameter<GibbsParameter>(p.first), p.second);
		}
	}

	return 0;
}

////////////////////////
// ConjugateParameter //
////////////////////////

double ConjugateParameter::_getNextSample(void) {
	return samplingFunction(_gibbs->getCurrentParameterValues());
}



/////////////////////////////
// ConstantVectorParameter //
/////////////////////////////

void ConstantVectorParameter::createElements(GibbsSampler* gibbs) {

	for (std::string elem : _elementNames) {
		gibbs->removeParameter(elem);
	}
	_elementNames.clear();

	for (unsigned int i = 0; i < fixedValues.size(); i++) {
		ConstantParameter cp;

		std::ostringstream oss;
		oss << i;

		cp.name = this->name + "[" + oss.str() + "]";
		cp.group = this->name + "_elem";
		cp.fixedValue = this->fixedValues[i];

		gibbs->addParameter(cp);

		_elementNames.push_back(cp.name);
	}
}

void ConstantVectorParameter::_parameterAddedToSampler(GibbsSampler* gs) {
	createElements(gs);
}

void ConstantVectorParameter::_parameterDeletedFromSampler(GibbsSampler* gs) {
	for (std::string elem : _elementNames) {
		gs->removeParameter(elem);
	}
}


////////////////////////
// VectorMH_Parameter //
////////////////////////
void VectorMH_Parameter::update(void) {
	//take samples
	std::vector<double> samples = __getNextSamples();

	//Regardless of whether there was an acceptance or rejection, samples contains the current values.
	for (unsigned int i = 0; i < _elementNames.size(); i++) {
		//store samples in elements
		UNSAFE_storeNewSample(_gibbs->getParameter<VectorElement>(_elementNames[i]), samples[i]);

		//update sample values. don't update self, update elements
		_gibbs->setCurrentParameterValue(_elementNames[i], samples[i]);
	}
}

void VectorMH_Parameter::_parameterAddedToSampler(GibbsSampler* gs) {
	createElements(gs);
}

void VectorMH_Parameter::_parameterDeletedFromSampler(GibbsSampler* gs) {
	for (std::string elem : _elementNames) {
		gs->removeParameter(elem);
	}
}

void VectorMH_Parameter::createElements(GibbsSampler* gibbs) {

	for (std::string elem : _elementNames) {
		gibbs->removeParameter(elem);
	}
	_elementNames.clear();

	for (unsigned int i = 0; i < startValues.size(); i++) {
		VectorElement ve;

		//ve.sourceName = this->name;
		//ve.index = i;

		std::ostringstream oss;
		oss << i;

		ve.name = this->name + "[" + oss.str() + "]";
		ve.group = this->name + "_elem";

		//Note: If an element with the same name already exists, you will be sad. Possibly check for existing?
		gibbs->addParameter(ve, startValues[i]);

		_elementNames.push_back(ve.name);
	}
}


std::vector<double> VectorMH_Parameter::__getNextSamples(void) {
	ParameterList param = _gibbs->getCurrentParameterValues();

	std::vector<double> currentValues(_elementNames.size());
	for (unsigned int i = 0; i < currentValues.size(); i++) {
		currentValues[i] = param[_elementNames[i]];
	}

	std::vector<double> candidateValues = this->deviateFunction(currentValues);

	for (unsigned int i = 0; i < candidateValues.size(); i++) {
		//if any of the candidates are outside of the range of possible values, reject all
		if (candidateValues[i] < ranges[i].first || candidateValues[i] > ranges[i].second) {
			acceptanceTracker.parameterRejected();
			return currentValues;
		}
	}

	//I probably don't need to do this, because the current values should already be in, but safety...
	for (unsigned int i = 0; i < currentValues.size(); i++) {
		param[_elementNames[i]] = currentValues[i];
	}
	_gibbs->updateDependentParameters(&param); // Probably don't need to here

	double currentLL = this->llFunction(currentValues, param); //Calculate current log likelihood


	//Update values for candidate values
	for (unsigned int i = 0; i < candidateValues.size(); i++) {
		param[_elementNames[i]] = candidateValues[i];
	}
	_gibbs->updateDependentParameters(&param);

	double candidateLL = this->llFunction(candidateValues, param); //Calculate candidate log likelihood


	double likelihoodRatio = exp(candidateLL - currentLL);

	std::vector<double> newValues;
	double canonicalSample = GibbsSampler::canonical(_gibbs->getGenerator());

	//if the likelihood ratio is greater than 1, the sample is always less than it, so accept
	if (canonicalSample < likelihoodRatio) {
		newValues = candidateValues;
		acceptanceTracker.parameterAccepted();
	} else {
		newValues = currentValues;
		acceptanceTracker.parameterRejected();
	}

	return newValues;
}

const std::vector<std::string>& VectorMH_Parameter::getElementNames(void) const {
	return _elementNames;
}

std::vector<double> VectorMH_Parameter::getCurrentValue(void) const {
	std::vector<double> values(_elementNames.size());
	ParameterList param = _gibbs->getParameterList(_elementNames);

	for (unsigned int i = 0; i < _elementNames.size(); i++) {
		values[i] = param.at(_elementNames[i]);
	}
	return values;
}

std::vector<double> VectorMH_Parameter::getIterationValue(unsigned int iteration) const {

	std::vector<GibbsParameter*> param = _gibbs->getParameters(_elementNames);

	std::vector<double> values(_elementNames.size());

	for (unsigned int i = 0; i < _elementNames.size(); i++) {
		values[i] = param[i]->getSample(iteration);
	}
	return values;
}

