#pragma once

#include "Compilation.h"
#include "UtilityFunctions.h"
#include "GibbsParameters.h"

class GibbsSampler;

struct SectionTracker {

	// Add section markers
	void sectionStart(std::string sectionName, std::function<void(std::string, bool)> callback = nullptr);
	void sectionEnd(std::string sectionName);

	// Clear section markers, effectively disabling it
	void clear(void);
	bool hasSections(void);

	void trackSectionTiming(bool track);
	double getTotalSecondsInSection(std::string sectionName);
	std::string getTimingResultsString(void);

private:

	struct SectionInfo {
		std::string name;
		std::function<void(std::string, bool)> callback;

		size_t startIndex;
		size_t endIndex;

		std::chrono::microseconds totalTimeInSection;
		std::chrono::time_point<std::chrono::steady_clock> sectionStartTime;
	};

	friend class GibbsSampler;

	SectionTracker(GibbsSampler* parent);

	void _startRun(void);
	void _update(size_t parameterIndex);

	GibbsSampler* _sampler;
	bool _trackTiming;

	std::vector<SectionInfo> _sections;
	std::vector<size_t> _sectionOrder;

};


class GibbsSampler {
public:

	SectionTracker sectionTracker;

	size_t iterationsPerStatusUpdate;

	//Takes a pointer to this and the current iteration. Returns true to continue, false to stop.
	std::function<bool(GibbsSampler*, size_t)> iterationCompleteCallback;

	GibbsSampler(void);

	void run(size_t samplesToCollect, bool clearExistingSamples = true);

	void clear(void);

	void createConstantParameter(std::string name, std::string group, double value, std::string replaceName = "", bool removeAndAdd = true);
	void setParameterGroupToConstantValue(std::string group, double value);
	void setParameterGroupDependency(std::string group, std::string sourceParameter, bool removeAndAdd = true);

	size_t getIndex(const std::string& parameterName) const;
	bool hasParameter(const std::string& parameterName) const;

	std::vector<GibbsParameter*> getParameters(void);
	std::vector<GibbsParameter*> getParameters(const std::vector<std::string>& param);

	std::vector<size_t> namesToIndices(std::vector<std::string> names) const;

	std::vector<double> getCurrentGroupValues(std::string group) const;


	void setCurrentParameterValue(std::string p, double v);
	const ParameterList& getCurrentParameterValues(void);
	ParameterList getParameterList(std::vector<std::string> parameterNames);
	ParameterList getParameterList(std::vector<size_t> parameterIndices);
	ParameterList getIterationParameterValues(size_t iteration);

	void updateDependentParameters(ParameterList* param) const;


	std::mt19937_64& getGenerator(void) {
		return _generator;
	}

	// Used internally, but users can use it as well
	static std::uniform_real_distribution<double> canonical;

	
	template <typename T>
	void addParameter(T param, double startValue = 0) {

		if (this->hasParameter(param.name)) {
			GS_COUT << "GibbsSampler::addParameter(): Warning: Attempt to add a parameter with a duplicate name. The name is: " << param.name << std::endl;
		}

		T* pp = new T;
		*pp = param;

		pp->_storeNewSample(startValue);
		pp->_gibbs = this;

		_paramVector.push_back(pp);

		size_t index = _paramVector.size() - 1;

		_groupToIndices[pp->group].push_back(index);

		_nameToIndex[param.name] = index;

		_paramVector[index]->_parameterAddedToSampler(this);
	}

	//replaceName = name of parameter to be replaced
	//This puts the new parameter in the place of the old parameter so it has the same update order
	template <typename T>
	void replaceParameter(std::string replaceName, T param, double startValue = 0) {

		if (_nameToIndex.find(replaceName) == _nameToIndex.end()) {
			GS_COUT << "Unable to replace parameter \"" << replaceName << "\" because it was not found." << std::endl;
			return;
		}

		//Can't you just use _remakeIndices() instead of all of this shit?
		//Delete the old name to index and make the new name point to that index
		int index = _nameToIndex[replaceName];
		_nameToIndex.erase(replaceName);
		_nameToIndex[param.name] = index;

		//Remove the old parameter from the _groupToIndices map
		std::string oldGroup = _paramVector[index]->group;
		std::vector<size_t>& gi = _groupToIndices[oldGroup];
		auto it = std::find(gi.begin(), gi.end(), index);
		if (it != gi.end()) {
			gi.erase(it);
		}

		//Add the new parameter to its _groupToIndices map
		_groupToIndices[param.group].push_back(index);


		//Add the new parameter
		T* pp = new T; //Allocate space for the new parameter
		*pp = param; //Copy the new parameter to that space
		pp->_storeNewSample(startValue);
		pp->_gibbs = this;

		_paramVector[index]->_parameterDeletedFromSampler(this); //Clean up old parameter
		delete _paramVector[index]; //Delete memory that was allocated for the old parameter

		_paramVector[index] = pp; //Store pointer to parameter copy
		_paramVector[index]->_parameterAddedToSampler(this);
	}

	bool removeParameter(std::string removeName) {
		if (_nameToIndex.find(removeName) == _nameToIndex.end()) {
			return false;
		}

		int index = _nameToIndex[removeName];

		_paramVector[index]->_parameterDeletedFromSampler(this);

		delete _paramVector[index]; //delete allocated memory
		_paramVector.erase(_paramVector.begin() + index);

		_remakeIndices();
		return true;
	}

	template <typename T>
	T* getParameter(std::string name) {
		return getParameter<T>(_nameToIndex.at(name));
	}

	template <typename T>
	const T* getParameter(std::string name) const {
		return getParameter<T>(_nameToIndex.at(name));
	}

	template <typename T>
	T* getParameter(size_t index) {
		GibbsParameter* p = _paramVector.at(index);
		return dynamic_cast<T*>(p);
	}

	template <typename T>
	const T* getParameter(size_t index) const {
		GibbsParameter* p = _paramVector.at(index);
		return dynamic_cast<const T*>(p);
	}

	template <typename T>
	std::vector<T*> getParameterGroup(std::string groupName) {
		std::vector<T*> p;
		for (GibbsParameter* par : _paramVector) {
			if (par->group == groupName) {
				p.push_back(dynamic_cast<T*>(par));
			}
		}
		return p;
	}


#ifdef COMPILING_WITH_CX
	void outputPosteriors(std::string outputDirectory);
	CX_DataFrame getPosteriors(void);

	CX_DataFrame getAcceptanceRates(void);
	void outputAcceptanceRates(std::string filename);
#endif
#ifdef COMPILING_WITH_RCPP
	Rcpp::DataFrame getAcceptanceRates(void);

	Rcpp::List getPosteriors(void);
#endif


private:

	std::map<std::string, std::vector<size_t>> _groupToIndices;

	std::map<std::string, size_t> _nameToIndex;
	std::vector<GibbsParameter*> _paramVector;

	ParameterList _currentValues;

	size_t _storedIterations;

	std::mt19937_64 _generator;
	
	void _remakeIndices(void);

	void _makeDependentParameterList(void);
	std::vector<size_t> _continuouslyUpdatedDependentParameters;
	std::vector<size_t> _dependentParameters;
	std::vector<size_t> _independentParameters;

};


