#include "GibbsSampler.h"


GibbsSampler::GibbsSampler(void) :
	sectionTracker(this),
	iterationsPerStatusUpdate(100),
	iterationCompleteCallback(nullptr),
	_storedIterations(0)
{
	std::random_device rd;
	_generator.seed(rd());
}


void GibbsSampler::createConstantParameter(std::string name, std::string group, double value, std::string replaceName, bool removeAndAdd) {
	DependentParameter dp;
	dp.name = name;
	dp.group = group;

	dp.sourceParameter = dp.name;

	if (replaceName == "") {
		this->addParameter(dp, value);
	} else {
		if (removeAndAdd) {
			this->removeParameter(replaceName);
			this->addParameter(dp, value);
		} else {
			this->replaceParameter(replaceName, dp, value);
		}
	}
}

ParameterList GibbsSampler::getParameterList(std::vector<std::string> parameterNames) {
	return getParameterList(namesToIndices(parameterNames));
}

ParameterList GibbsSampler::getParameterList(std::vector<size_t> parameterIndices) {
	ParameterList rval;

	for (size_t i = 0; i < parameterIndices.size(); i++) {
		rval[_paramVector.at(i)->name] = _paramVector.at(i)->value();
	}

	return rval;
}

const ParameterList& GibbsSampler::getCurrentParameterValues(void) {
	return _currentValues;
}

void GibbsSampler::setCurrentParameterValue(std::string p, double v) {
	_currentValues.at(p) = v;
}

ParameterList GibbsSampler::getIterationParameterValues(size_t iteration) {

	ParameterList rval;

	if (iteration >= _storedIterations) {
		return rval;
	}

	for (size_t j = 0; j < _paramVector.size(); j++) {
		GibbsParameter* p = _paramVector[j];

		switch (p->getDimension()) {
		case ParameterDimension::SCALAR:
			rval[p->name] = p->getSample(iteration);
			break;
		case ParameterDimension::VECTOR:
		{
			std::map<std::string, double> values = dynamic_cast<BaseVectorParameter*>(p)->getIterationValueMap(iteration);
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
void GibbsSampler::run(size_t samplesToCollect, bool clearExistingSamples) {

	_makeDependentParameterList();

	for (size_t j = 0; j < _paramVector.size(); j++) {
		GibbsParameter* p = _paramVector[j];

		//By using getSamples(), it deals with the fact that constant parameters do not have a value stored in _samples (but what about VectorElement?).
		size_t samples = p->getSamples().size();

		if (samples == 0) {
			GS_COUT << "Warning: Parameter " << p->name << " has no starting value. It has had 0 added as the starting value." << std::endl;
			p->_samples.push_back(0);
		}
	}

	// At least one iteration is stored (the starting values)
	_storedIterations = std::max<size_t>(_storedIterations, 1);

	if (clearExistingSamples) {
		for (size_t j = 0; j < _paramVector.size(); j++) {
			//clear all but the starting value
			_paramVector[j]->_samples.resize(1); //TODO: Bad use of private vars

			if (_paramVector[j]->getAcceptanceTracker() != nullptr) {
				_paramVector[j]->getAcceptanceTracker()->reset();
			}
		}

		_storedIterations = 1; //Just the starting values
	}

	
	for (size_t j = 0; j < _paramVector.size(); j++) {
		// Default initialization, to be overridden.
		_currentValues[_paramVector[j]->name] = 0;

		// Make sure the parameters are pointing to this.
		_paramVector[j]->_gibbs = this;
	}

	// Update the current values for the parameters
	// Do not use setCurrentParameterValue because that function assumes that parameters already exist.
	for (size_t j = 0; j < _independentParameters.size(); j++) {
		size_t ind = _independentParameters[j];
		_currentValues[_paramVector[ind]->name] = _paramVector[ind]->value();
	}
	for (size_t j = 0; j < _dependentParameters.size(); j++) {
		size_t ind = _dependentParameters[j];
		_currentValues[_paramVector[ind]->name] = _paramVector[ind]->value();
	}
	

	auto runStartTime = std::chrono::high_resolution_clock::now();
	auto lapStartTime = runStartTime;

	sectionTracker.trackSectionTiming(true); // TODO: Make configurable
	sectionTracker._startRun();

	for (size_t i = 0; i < samplesToCollect; i++) {
		
		//Update the parameters
		for (size_t j = 0; j < _paramVector.size(); j++) {
			
			// Update sections before parameters
			sectionTracker._update(j);

			// Update the parameter and dependent parameters
			_paramVector[j]->update();
			updateDependentParameters(&_currentValues);

		}

		// Update one last time for the past-the-end index
		sectionTracker._update(_paramVector.size());

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

			auto currentTime = std::chrono::high_resolution_clock::now();
			long long lapDurationMs = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lapStartTime).count();
			long long msPerIteration = (double)lapDurationMs / iterationsPerStatusUpdate;

			size_t remainingIterations = samplesToCollect - i;
			long long minRemaining = remainingIterations * msPerIteration / (1000.0 * 60.0);

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

#if COMPILING_WITH_RCPP
			//Check for interrupts whenever a status update happens.
			Rcpp::checkUserInterrupt();
#endif
		}
	}
	// run complete

#if COMPILING_WITH_RCPP
	GS_COUT << "\r100% complete.                                                      " << std::endl;
#endif

	if (sectionTracker.hasSections()) {
		GS_COUT << sectionTracker.getTimingResultsString() << std::endl;
	}
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
	for (size_t i = 0; i < values.size(); i++) {
		ss << values[i] << std::endl;
	}
	Util::writeToFile(filename, ss.str(), false);
}

void GibbsSampler::outputPosteriors(std::string outputDirectory) {

	for (GibbsParameter* p : _paramVector) {
		GS_COUT << "Outputting posterior for " << p->name << std::endl;
		outputDataVector(outputDirectory + "/" + p->name + ".txt", p->name, p->getSamples());
	}

}

//This function is probably really slow.
CX_DataFrame GibbsSampler::getPosteriors(void) {

	CX_DataFrame post;

	for (GibbsParameter* p : _paramVector) {
		vector<double> samp = p->getSamples();
		for (size_t i = 0; i < samp.size(); i++) {
			post(i, p->name) = samp[i];
		}
	}

	return post;
}

#endif

#ifdef COMPILING_WITH_RCPP

Rcpp::DataFrame GibbsSampler::getAcceptanceRates(void) {

	std::vector<GibbsParameter*> param = this->getParameters();

	Rcpp::StringVector names(param.size());
	Rcpp::StringVector groups(param.size());
	Rcpp::NumericVector acceptanceRates(param.size());

	for (size_t i = 0; i < param.size(); i++) {

		names[i] = param[i]->name;
		groups[i] = param[i]->group;

		if (param[i]->getAcceptanceTracker() != nullptr) {
			acceptanceRates[i] = param[i]->getAcceptanceTracker()->acceptanceRate();
		} else {
			acceptanceRates[i] = NA_REAL; // I don't know where this symbol comes from, but it magically appears.
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

	for (size_t i = 0; i < param.size(); i++) {
	
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

	for (size_t i = 0; i < _paramVector.size(); i++) {
		_nameToIndex[_paramVector[i]->name] = i;
		_groupToIndices[_paramVector[i]->group].push_back(i);
	}
}

void GibbsSampler::_makeDependentParameterList(void) {
	_continuouslyUpdatedDependentParameters.clear();
	_dependentParameters.clear();
	_independentParameters.clear();

	for (size_t i = 0; i < _paramVector.size(); i++) {
		GenericDependentParameter* p = (GenericDependentParameter*)_paramVector[i];
		if (p->hasDependency()) {
			if (p->updateContinuously) {
				_continuouslyUpdatedDependentParameters.push_back(i);
			}
			_dependentParameters.push_back(i);
		} else {
			_independentParameters.push_back(i);
		}
	}
}

/*
Recalculates values for all continuously updated dependent parameters, in their insertion order.
*/
void GibbsSampler::updateDependentParameters(ParameterList* param) const {
	for (size_t i : _continuouslyUpdatedDependentParameters) {
		GenericDependentParameter* p = (GenericDependentParameter*)_paramVector[i];
		(*param)[p->name] = p->evaluate(*param);
	}
}


/* Replaces all of the parameters in the given parameter group with a ConstantParameter with the given value.
*/
void GibbsSampler::setParameterGroupToConstantValue(std::string group, double value) {
	std::vector<GibbsParameter*> params = this->getParameterGroup<GibbsParameter>(group);
	for (size_t i = 0; i < params.size(); i++) {
		ConstantParameter cp(value, params[i]->name, params[i]->group);
		this->replaceParameter(params[i]->name, cp, value);
	}
}

void GibbsSampler::setParameterGroupDependency(std::string group, std::string sourceParameter, bool removeAndAdd) {
	std::vector<GibbsParameter*> params = this->getParameterGroup<GibbsParameter>(group);
	for (size_t i = 0; i < params.size(); i++) {
		DependentParameter dp(params[i]->name, params[i]->group, sourceParameter);

		if (removeAndAdd) {
			this->removeParameter(params[i]->name);
			this->addParameter(dp);
		} else {
			this->replaceParameter(params[i]->name, dp);
		}
	}
}

std::vector<size_t> GibbsSampler::namesToIndices(std::vector<std::string> names) const {
	std::vector<size_t> indices(names.size());
	for (size_t i = 0; i < indices.size(); i++) {
		indices[i] = _nameToIndex.at(names[i]);
	}
	return indices;
}

std::vector<double> GibbsSampler::getCurrentGroupValues(std::string group) const {
	std::vector<size_t> groupIndices = _groupToIndices.at(group);
	std::vector<double> rval(groupIndices.size());
	for (size_t i = 0; i < groupIndices.size(); i++) {
		rval[i] = _paramVector[groupIndices[i]]->value();
	}
	return rval;
}

/* Returns -1 if the parameter is not found. */
size_t GibbsSampler::getIndex(const std::string& parameterName) const {
	if (_nameToIndex.find(parameterName) == _nameToIndex.end()) {
		return (size_t)-1;
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

	for (size_t i = 0; i < rval.size(); i++) {

		size_t index = this->getIndex(param[i]);

		if (index != (size_t)-1) {
			rval[i] = _paramVector[index];
		} else {
			rval[i] = nullptr;
		}
	}

	return rval;
}

void GibbsSampler::clear(void) {
	for (size_t i = 0; i < _paramVector.size(); i++) {
		delete _paramVector[i];
	}
	_paramVector.clear();
	_nameToIndex.clear();

	_storedIterations = 0;

	sectionTracker.clear();

}

////////////////////
// SectionTracker //
////////////////////

SectionTracker::SectionTracker(GibbsSampler* gibbs) :
	_sampler(gibbs),
	_trackTiming(false)
{
	this->clear();
}


void SectionTracker::sectionStart(std::string sectionName, std::function<void(std::string, bool)> callback) {

	SectionInfo info;
	info.name = sectionName;
	info.callback = callback;

	// Section tracker update happens before parameter update
	size_t nextParamIndex = _sampler->getParameters().size();
	info.startIndex = nextParamIndex;

	_sections.push_back(info);

	//_sectionOrder.push_back(_sections.size() - 1);

}

void SectionTracker::sectionEnd(std::string sectionName) {

	size_t i = 0; // section index of section start

	for (; i < _sections.size(); i++) {
		if (_sections[i].name == sectionName) {
			break;
		}
	}
	if (i >= _sections.size()) {
		GS_COUT << "Error adding section end: Section name: \"" << sectionName << "\" not found.";
		return;
	}

	size_t nextParamIndex = _sampler->getParameters().size();
	_sections[i].endIndex = nextParamIndex;

	//_sectionOrder.push_back(i);

}

void SectionTracker::_startRun(void) {

	for (SectionInfo& sec : _sections) {
		sec.totalTimeInSection = std::chrono::microseconds(0);
	}

}

// Sections can start and end on the same index,
// have a start but no end, etc. This update is inefficient, but simple
void SectionTracker::_update(size_t index) {

	for (SectionInfo& sec : _sections) {

		if (sec.startIndex == index) {

			if (sec.callback) {
				sec.callback(sec.name, true);
			}

			if (_trackTiming) {
				sec.sectionStartTime = std::chrono::steady_clock::now();
			}
		}

		if (sec.endIndex == index) {

			if (_trackTiming) {
				auto endTime = std::chrono::steady_clock::now();
				std::chrono::microseconds newMicros = std::chrono::duration_cast<std::chrono::microseconds>(endTime - sec.sectionStartTime);
				sec.totalTimeInSection += newMicros;
			}

			if (sec.callback) {
				sec.callback(sec.name, false);
			}

		}
	}
}

void SectionTracker::clear(void) {

	_sections.clear();
	//_sectionOrder.clear();
}

bool SectionTracker::hasSections(void) {
	return _sections.size() > 0;
}

void SectionTracker::trackSectionTiming(bool track) {
	_trackTiming = track;
}

// in seconds
double SectionTracker::getTotalSecondsInSection(std::string sectionName) {

	for (SectionInfo& sec : _sections) {
		if (sec.name == sectionName) {
			return sec.totalTimeInSection.count() / 1000000.0;
		}
	}

	return 0;
}

std::string SectionTracker::getTimingResultsString(void) {

	std::vector<double> secs(_sections.size());
	double totalSecs = 0;

	for (size_t i = 0; i < _sections.size(); i++) {

		secs[i] = _sections[i].totalTimeInSection.count() / 1000000.0;
		totalSecs += secs[i];

	}

	std::ostringstream oss;

	oss << "Section Timing Results:" << std::endl;

	for (size_t i = 0; i < _sections.size(); i++) {

		oss << std::setw(4) << (int)(secs[i] / totalSecs * 100.0) << "  ";

		oss << _sections[i].name << " (" << secs[i] << " sec)";

		oss << std::endl;

	}

	return oss.str();
}