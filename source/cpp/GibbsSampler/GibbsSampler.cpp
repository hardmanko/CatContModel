#include "GibbsSampler.h"

std::uniform_real_distribution<double> GibbsSampler::canonical = std::uniform_real_distribution<double>(0, 1);

//This is here because it uses the GibbsSampler class
void GibbsParameter::_updateCurrentValue(double s) {
	_gibbs->getCurrentParameterValues().at(this->name) = s;
}

//////////////////
// MH_Parameter //
//////////////////
double MH_Parameter::_getNextSample(void) {

	ParameterList& param = _gibbs->getCurrentParameterValues();

	double currentValue = this->value();
	double candidate = this->deviateFunction(currentValue);

	//if the candidate is outside of the range of possible values, reject it
	if (candidate < range.lower || candidate > range.upper) {
		return currentValue;
	}
	
	param[this->name] = currentValue;
	double currentLL = this->llFunction(currentValue, param);

	param[this->name] = candidate;
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
	double llCandidate = llFunction(candidateValues, paramCopy);


	double likelihoodRatio = exp(llCandidate - llCurrent);

	//double newValue;
	double canonicalSample = GibbsSampler::canonical(_gibbs->getGenerator());

	ParameterList& paramToUpdate = _gibbs->getCurrentParameterValues();

	//if the likelihood ratio is greater than 1, the sample is always less than it, so accept
	if (canonicalSample < likelihoodRatio) {

		for (const auto& p : candidateValues) {
			UNSAFE_storeNewSample(_gibbs->getParameter<GibbsParameter>(p.first), p.second);

			paramToUpdate.at(p.first) = p.second; //if accepting, change the current values of the parameters stored by the sampler
			//if not accepting, nothing needs to be done to the current values as they have not been changed
		}

		acceptanceTracker.parameterAccepted();
	} else {
		for (const auto& p : currentValues) {
			UNSAFE_storeNewSample(_gibbs->getParameter<GibbsParameter>(p.first), p.second);
		}
		acceptanceTracker.parameterRejected();
	}

	return 0;
}



double ConjugateParameter::_getNextSample(void) {
	return samplingFunction(_gibbs->getCurrentParameterValues());
}

double DependentParameter::_getNextSample(void) {
	const ParameterList& pl = _gibbs->getCurrentParameterValues();
	if (sourceParameter != "") {
		return pl.at(sourceParameter);
	} else if (samplingFunction != nullptr) {
		return samplingFunction(pl);
	}
	throw(std::runtime_error("DependentParameter not configured with either a parameter name or a sampling function!"));
	return 0;
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

	ParameterList& param = _gibbs->getCurrentParameterValues();


	for (unsigned int i = 0; i < _elementNames.size(); i++) {
		//store samples in elements
		UNSAFE_storeNewSample(_gibbs->getParameter<VectorElement>(_elementNames[i]), samples[i]);

		//update sample values. don't update self, update elements
		param[_elementNames[i]] = samples[i];
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
	ParameterList& param = _gibbs->getCurrentParameterValues();

	std::vector<double> currentValues(_elementNames.size());
	for (unsigned int i = 0; i < currentValues.size(); i++) {
		currentValues[i] = param[_elementNames[i]];
	}

	std::vector<double> candidateValues = this->deviateFunction(currentValues);

	for (unsigned int i = 0; i < candidateValues.size(); i++) {
		//if any of the candidates are outside of the range of possible values, reject all
		if (candidateValues[i] < ranges[i].first || candidateValues[i] > ranges[i].second) {
			return currentValues;
		}
	}

	//I probably don't need to do this, because the current values should already be in, but safety...
	for (unsigned int i = 0; i < currentValues.size(); i++) {
		param[_elementNames[i]] = currentValues[i];
	}

	//Calculate current log likelihood
	double currentLL = this->llFunction(currentValues, param);

	//Update values for candidate values
	for (unsigned int i = 0; i < candidateValues.size(); i++) {
		param[_elementNames[i]] = candidateValues[i];
	}

	//calculate candidate log likelihood
	double candidateLL = this->llFunction(candidateValues, param);


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
	ParameterList& param = _gibbs->getCurrentParameterValues();
	for (unsigned int i = 0; i < _elementNames.size(); i++) {
		values[i] = param[_elementNames[i]];
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


//////////////////
// GibbsSampler //
//////////////////
ParameterList GibbsSampler::getParameterList(std::vector<std::string> parameterNames) {
	return getParameterList(namesToIndices(parameterNames));
}

ParameterList GibbsSampler::getParameterList(std::vector<unsigned int> parameterIndices) {
	ParameterList rval;

	for (unsigned int i = 0; i < parameterIndices.size(); i++) {
		rval[_paramVector.at(i)->name] = _paramVector.at(i)->value();
	}

	return rval;
}

ParameterList& GibbsSampler::getCurrentParameterValues(void) {
	return _currentValues;
}

ParameterList GibbsSampler::getIterationParameterValues(unsigned int iteration) {

	ParameterList rval;

	if (iteration >= _storedIterations) {
		return rval;
	}

	for (unsigned int j = 0; j < _paramVector.size(); j++) {
		GibbsParameter* p = _paramVector[j];

		switch (p->dimension) {
		case ParameterDimension::SCALAR:
			rval[p->name] = p->getSample(iteration);
			break;
		case ParameterDimension::VECTOR:
		{
			std::map<std::string, double> values = dynamic_cast<VectorParameter*>(p)->getIterationValueMap(iteration);
			for (auto it = values.begin(); it != values.end(); it++) {
				rval[it->first] = it->second;
			}
		}
			break;
		case ParameterDimension::OTHER:
		case ParameterDimension::ZERO:
			break;
		}

	}

	return rval;

}

/*

clearExistingSamples clears all samples except for the starting value. It also resets MH acceptance rates.
*/
void GibbsSampler::run(unsigned int samplesToCollect, bool clearExistingSamples) {

	for (unsigned int j = 0; j < _paramVector.size(); j++) {
		GibbsParameter* p = _paramVector[j];

		//By using getSamples(), it deals with the fact that constant parameters do not have a value stored in _samples (but what about VectorElement?).
		unsigned int samples = p->getSamples().size();

		if (samples == 0) {
			GS_COUT << "Warning: Parameter " << p->name << " has no starting value. It has had 0 added as the starting value." << std::endl;
			p->_samples.push_back(0);
		}
	}

	//At least one iteration is stored (the starting values)
	_storedIterations = std::max<unsigned int>(_storedIterations, 1);

	if (clearExistingSamples) {
		for (unsigned int j = 0; j < _paramVector.size(); j++) {
			//clear all but the starting value
			_paramVector[j]->_samples.resize(1);

			if (_paramVector[j]->getAcceptanceTracker() != nullptr) {
				_paramVector[j]->getAcceptanceTracker()->reset();
			}
		}

		_storedIterations = 1; //Just the starting values
	}

	
	for (unsigned int j = 0; j < _paramVector.size(); j++) {
		//Update the current values from the parameters
		_currentValues[_paramVector[j]->name] = _paramVector[j]->value();

		//make sure the parameters are pointing to this.
		_paramVector[j]->_gibbs = this;
	}


	auto runStartTime = std::chrono::high_resolution_clock::now();
	auto lapStartTime = runStartTime;

	for (unsigned int i = 0; i < samplesToCollect; i++) {
		
		//Update the parameters
		for (unsigned int j = 0; j < _paramVector.size(); j++) {
			_paramVector[j]->update();
		}

		_storedIterations++;

		//Do the iteration status update function call
		if (iterationCompleteCallback) {
			bool continue_ = iterationCompleteCallback(this, i);
			if (!continue_) {
				GS_COUT << "Stopped by iteration complete callback." << std::endl;
				return;
			}
		}

		//Do status update
		if ((i % iterationsPerStatusUpdate) == 0) {

			//CX_Millis elapsedTime = Clock.now() - lapStartTime;
			auto currentTime = std::chrono::high_resolution_clock::now();
			long long lapDurationMs = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lapStartTime).count();
			long long msPerIteration = (double)lapDurationMs / iterationsPerStatusUpdate;

			long long totalDurationSec = std::chrono::duration_cast<std::chrono::seconds>(currentTime - runStartTime).count();
			double proportionComplete = ((double)i + 1) / samplesToCollect;
			double proportionRemaining = 1 - proportionComplete;
			double secondsRemaining = proportionRemaining / proportionComplete * totalDurationSec;
			long long minRemaining = secondsRemaining / 60.0;

#if COMPILING_WITH_RCPP
			GS_COUT << "\r"; // Does not print new line, just overwrites old line.
#endif
			GS_COUT << 100.0f * (float)i / samplesToCollect << "% complete. " <<
				msPerIteration << " ms per iteration. " <<
				minRemaining << " min remaining.";
#if COMPILING_WITH_RCPP
			GS_COUT << "            "; // Extra spaces cover for different length lines.
#else
			GS_COUT << std::endl;
#endif

			lapStartTime = std::chrono::high_resolution_clock::now();
		}
	}

#if COMPILING_WITH_RCPP
	GS_COUT << "\r100% complete.                                                      " << std::endl;
#endif
}



#ifdef COMPILING_WITH_CX

void GibbsSampler::outputAcceptanceRates(std::string filename) {
	CX_DataFrame accept = getAcceptanceRates();
	accept.printToFile(filename);
}

CX_DataFrame GibbsSampler::getAcceptanceRates(void) {
	CX_DataFrame accept;

	std::vector<GibbsParameter*> param = this->getParameters();
	for (auto p : param) {
		CX_DataFrameRow row;
		row["name"] = p->name;
		row["group"] = p->group;
		if (p->getAcceptanceTracker() != nullptr) {
			row["acceptanceRate"] = p->getAcceptanceTracker()->acceptanceRate();
		} else {
			row["acceptanceRate"] = "NA";
		}

		accept.appendRow(row);
	}

	return accept;
}

void outputDataVector(std::string filename, std::string header, std::vector<double> values) {
	std::ostringstream ss;
	ss << std::scientific << std::setprecision(20);
	ss << header << std::endl;
	for (unsigned int i = 0; i < values.size(); i++) {
		ss << values[i] << std::endl;
	}
	Util::writeToFile(filename, ss.str(), false);
}

void GibbsSampler::outputPosteriors(std::string outputDirectory) {

	std::vector<GibbsParameter*> param = this->getParameters();
	for (GibbsParameter* p : param) {
		GS_COUT << "Outputting posterior for " << p->name << std::endl;
		outputDataVector(outputDirectory + "/" + p->name + ".txt", p->name, p->getSamples());
	}

}

#endif

#ifdef COMPILING_WITH_RCPP

Rcpp::DataFrame GibbsSampler::getAcceptanceRates(void) {

	std::vector<GibbsParameter*> param = this->getParameters();

	Rcpp::StringVector names(param.size());
	Rcpp::StringVector groups(param.size());
	Rcpp::NumericVector acceptanceRates(param.size());

	for (unsigned int i = 0; i < param.size(); i++) {

		names[i] = param[i]->name;
		groups[i] = param[i]->group;

		if (param[i]->getAcceptanceTracker() != nullptr) {
			acceptanceRates[i] = param[i]->getAcceptanceTracker()->acceptanceRate();
		} else {
			acceptanceRates[i] = NA_REAL; //This is really fucked up. I don't know where this symbol comes from, but it magically appears.
		}
	}

	Rcpp::DataFrame accept = Rcpp::DataFrame::create(Rcpp::Named("name") = names, 
		Rcpp::Named("group") = groups, 
		Rcpp::Named("acceptanceRate") = acceptanceRates);

	return accept;
}

Rcpp::List GibbsSampler::getPosteriors(void) {

	std::vector<GibbsParameter*> param = this->getParameters();

	Rcpp::List rval(param.size());
	Rcpp::CharacterVector names(param.size());

	for (unsigned int i = 0; i < param.size(); i++) {
	
		const std::vector<double> s = param[i]->getSamples();
		Rcpp::NumericVector v(s.begin(), s.end());
		rval[i] = v;

		names[i] = param[i]->name;

	}
	rval.attr("names") = names;

	return rval;
}

#endif


void GibbsSampler::_remakeIndices(void) {
	_groupToIndices.clear();
	_nameToIndex.clear();

	for (unsigned int i = 0; i < _paramVector.size(); i++) {
		_nameToIndex[_paramVector[i]->name] = i;
		_groupToIndices[_paramVector[i]->group].push_back(i);
	}
}

/* Replaces all of the parameters in the given parameter group with a ConstantParameter with the given value.
*/
void GibbsSampler::setParameterGroupToConstantValue(std::string group, double value) {
	std::vector<GibbsParameter*> params = this->getParameterGroup<GibbsParameter>(group);
	for (unsigned int i = 0; i < params.size(); i++) {
		std::string name = params[i]->name;
		this->replaceParameter(name, ConstantParameter(value, name, params[i]->group), value);
	}
}

std::vector<unsigned int> GibbsSampler::namesToIndices(std::vector<std::string> names) const {
	std::vector<unsigned int> indices(names.size());
	for (unsigned int i = 0; i < indices.size(); i++) {
		indices[i] = _nameToIndex.at(names[i]);
	}
	return indices;
}

std::vector<double> GibbsSampler::getCurrentGroupValues(std::string group) const {
	std::vector<unsigned int> groupIndices = _groupToIndices.at(group);
	std::vector<double> rval(groupIndices.size());
	for (unsigned int i = 0; i < groupIndices.size(); i++) {
		rval[i] = _paramVector[groupIndices[i]]->value();
	}
	return rval;
}

/* Returns -1 if the parameter is not found. */
unsigned int GibbsSampler::getIndex(const std::string& parameterName) const {
	if (_nameToIndex.find(parameterName) == _nameToIndex.end()) {
		return (unsigned int)-1;
	}

	return _nameToIndex.at(parameterName);
}

bool GibbsSampler::hasParameter(const std::string& parameterName) const {
	return _nameToIndex.find(parameterName) != _nameToIndex.end();
}

std::vector<GibbsParameter*> GibbsSampler::getParameters(void) {
	return _paramVector;
}

std::vector<GibbsParameter*> GibbsSampler::getParameters(const std::vector<std::string>& param) {

	std::vector<GibbsParameter*> rval(param.size());

	for (unsigned int i = 0; i < rval.size(); i++) {

		unsigned int index = this->getIndex(param[i]);

		if (index != (unsigned int)-1) {
			rval[i] = _paramVector[index];
		} else {
			rval[i] = nullptr;
		}
	}

	return rval;
}

void GibbsSampler::clear(void) {
	for (unsigned int i = 0; i < _paramVector.size(); i++) {
		delete _paramVector[i];
	}
	_paramVector.clear();
	_nameToIndex.clear();

	_storedIterations = 0;
}
