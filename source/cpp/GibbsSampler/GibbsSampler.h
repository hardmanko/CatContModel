#pragma once

#include <string>
#include <map>
#include <vector>
#include <functional>
#include <random>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <chrono>

#include "Compilation.h"

#if COMPILING_WITH_RCPP
#define GS_COUT Rcpp::Rcout
#endif
#if COMPILING_WITH_CX
#define GS_COUT std::cout
#endif


class GibbsSampler; //Forward declaration

typedef std::map<std::string, double> ParameterList; //A map string -> double is a conventient way to pass around parameter values

enum class ParameterDimension {
	ZERO,
	SCALAR,
	VECTOR,
	OTHER
};

class AcceptanceTracker {
public:

	AcceptanceTracker(void) :
		_acceptances(0),
		_rejections(0)
	{}

	void parameterAccepted(void) {
		_acceptances++;
	}

	void parameterRejected(void) {
		_rejections++;
	}

	unsigned int& acceptances(void) {
		return _acceptances;
	}

	unsigned int& rejections(void) {
		return _rejections;
	}

	double acceptanceRate(void) const {
		return ((double)_acceptances) / (_acceptances + _rejections);
	}

	void reset(void) {
		_acceptances = 0;
		_rejections = 0;
	}

private:
	unsigned int _acceptances;
	unsigned int _rejections;
};



class GibbsParameter {
public:

	GibbsParameter(void) :
		name("NULL"),
		group("NULL"),
		_gibbs(nullptr),
		_dimension(ParameterDimension::OTHER)
	{}

	virtual ~GibbsParameter(void) {
		//do nothing
	}

	std::string name;
	std::string group;

	
	ParameterDimension getDimension(void) {
		return _dimension;
	}

	virtual std::string type(void) const {
		return "GibbsParameter";
	}

	virtual void update(void) {
		double s = this->_getNextSample();

		this->_storeNewSample(s);

		this->_updateCurrentValue(s);
	}

	virtual double value(void) const {
		return _samples.back();
	}

	virtual double getSample(unsigned int iteration) const {
		//return this->getSamples.at(iteration);
		return _samples.at(iteration);
	}

	virtual std::vector<double>& getSamples(void) {
		return _samples;
	}

	virtual AcceptanceTracker* getAcceptanceTracker(void) {
		return nullptr;
	}

	virtual bool hasDependency(void) {
		return false;
	}

protected:

	virtual double _getNextSample(void) = 0; //This is the only pure virtual function

	virtual void _storeNewSample(double s) {
		_samples.push_back(s);
	}

	virtual void _updateCurrentValue(double s);

	virtual void _parameterAddedToSampler(GibbsSampler* gs) {
		return;
	}

	virtual void _parameterDeletedFromSampler(GibbsSampler* gs) {
		return;
	}


	friend class GibbsSampler;
	GibbsSampler* _gibbs;

	friend void UNSAFE_storeNewSample(GibbsParameter* gp, double s) {
		gp->_storeNewSample(s);
	}

	std::vector<double> _samples;
	ParameterDimension _dimension;
};



class VectorParameter : public GibbsParameter {
public:

	//pure virtual functions
	virtual void createElements(GibbsSampler* gs) = 0;

	virtual std::vector<double> getCurrentValue(void) const = 0;

	virtual std::vector<double> getIterationValue(unsigned int iteration) const = 0;



	virtual const std::vector<std::string>& getElementNames(void) const {
		return _elementNames;
	}


	virtual std::map<std::string, double> getCurrentValueMap(void) const {
		//map from element name to value
		const std::vector<std::string>& elem = getElementNames();
		std::vector<double> value = getCurrentValue();

		std::map<std::string, double> m;
		for (unsigned int i = 0; i < elem.size(); i++) {
			m[elem[i]] = value[i];
		}
		return m;
	}

	virtual std::map<std::string, double> getIterationValueMap(unsigned int iteration) const {
		//map from element name to value
		const std::vector<std::string>& elem = getElementNames();
		std::vector<double> value = getIterationValue(iteration);

		std::map<std::string, double> m;
		for (unsigned int i = 0; i < elem.size(); i++) {
			m[elem[i]] = value[i];
		}
		return m;
	}

protected:

	std::vector<std::string> _elementNames;

};



class ConjugateParameter : public GibbsParameter {
public:

	ConjugateParameter(void) {
		_dimension = ParameterDimension::SCALAR;
	}

	std::function<double(const ParameterList&)> samplingFunction;

	std::string type(void) const OVERRIDE {
		return "ConjugateParameter";
	}

private:

	double _getNextSample(void) OVERRIDE;

};

class GenericDependentParameter : public GibbsParameter {
public:

	GenericDependentParameter(void) {
		updateContinuously = true;
	}

	bool updateContinuously;

	bool hasDependency(void) OVERRIDE {
		return true;
	}

	virtual double evaluate(const ParameterList& param) const = 0;
};

class DependentParameter : public GenericDependentParameter {
public:

	DependentParameter(void) {
		_dimension = ParameterDimension::SCALAR;
	}

	std::string sourceParameter; //The name of the parameter from which to take its value

	std::string type(void) const OVERRIDE {
		return "DependentParameter";
	}

	double evaluate(const ParameterList& param) const OVERRIDE;
	double value(void) const OVERRIDE;

	double getSample(unsigned int iteration) const OVERRIDE;
	std::vector<double>& getSamples(void) OVERRIDE;

	virtual bool hasDependency(void) OVERRIDE {
		return true;
	}

protected:

	double _getNextSample(void) OVERRIDE;

};

class CalculatedParameter : public GenericDependentParameter {
public:
	CalculatedParameter(void) {
		_dimension = ParameterDimension::SCALAR;
	}

	//An arbitrary function that can do whatever (but probably calculate the parameter based on other parameters).
	std::function<double(const ParameterList&)> samplingFunction; 

	std::string type(void) const OVERRIDE {
		return "CalculatedParameter";
	}

	double evaluate(const ParameterList& param) const OVERRIDE;
	double value(void) const OVERRIDE;

	double getSample(unsigned int iteration) const OVERRIDE;
	std::vector<double>& getSamples(void) OVERRIDE;

	virtual bool hasDependency(void) OVERRIDE {
		return true;
	}

protected:

	double _getNextSample(void) OVERRIDE;
};

class ConstantParameter : public GibbsParameter {
public:

	ConstantParameter(void) :
		fixedValue(0)
	{
		_dimension = ParameterDimension::SCALAR;
	}

	ConstantParameter(double val, std::string name_, std::string group_) :
		fixedValue(val)
	{
		_dimension = ParameterDimension::SCALAR;

		this->name = name_;
		this->group = group_;
	};

	double fixedValue;

	double value(void) const OVERRIDE {
		return fixedValue;
	}

	std::string type(void) const OVERRIDE {
		return "ConstantParameter";
	}

	std::vector<double>& getSamples(void) OVERRIDE {
		_samples.clear();
		_samples.push_back(fixedValue);
		return _samples;
	}

	double getSample(unsigned int iteration) const OVERRIDE {
		return fixedValue;
	}

private:
	friend class GibbsSampler;

	double _getNextSample(void) OVERRIDE {
		return fixedValue;
	}

	void _storeNewSample(double s) OVERRIDE {
		return;
	}
};

//This class creates elements when it is added to the Gibbs Sampler, so all you have to do is
//call the convenience constructor.
class ConstantVectorParameter : public VectorParameter {
public:

	ConstantVectorParameter(const std::vector<double>& vals, std::string name_, std::string group_) {
		this->name = name_;
		this->group = group_;
		this->fixedValues = vals;
		_dimension = ParameterDimension::VECTOR;
	}

	std::vector<double> fixedValues;

	void update(void) OVERRIDE {
		return;
	}

	double value(void) const OVERRIDE {
		return std::numeric_limits<double>::quiet_NaN();
	}

	std::string type(void) const OVERRIDE {
		return "ConstantVectorParameter";
	}


	void createElements(GibbsSampler* gibbs) OVERRIDE;

	const std::vector<std::string>& getElementNames(void) const OVERRIDE {
		return _elementNames;
	}

	std::vector<double> getCurrentValue(void) const OVERRIDE {
		return fixedValues;
	}

	std::vector<double> getIterationValue(unsigned int iteration) const OVERRIDE {
		return fixedValues;
	}

private:

	double _getNextSample(void) OVERRIDE {
		return 0; // This is never called
	}

	void _parameterAddedToSampler(GibbsSampler* gs) OVERRIDE;

	void _parameterDeletedFromSampler(GibbsSampler* gs) OVERRIDE;

};

/*
The user of this class must set the following properies of this class:
* deviateFunction is a function that takes the value of the parameter for the previous iteration and returns a new candidate value for the parameter.
* llFunction takes the value of the parameter and a ParameterList of all of the other parameters in the 
model and returns the (potentially proportional) log likelihood for the parameter value, including the
prior on the parameter.
* The allowed range of the parameter may be optionally set with the range property. Any candidates outside of this range will be rejected.
*/
class MH_Parameter : public GibbsParameter {
public:

	MH_Parameter(void) {
		range.lower = -std::numeric_limits<double>::infinity();
		range.upper = std::numeric_limits<double>::infinity();

		_dimension = ParameterDimension::SCALAR;
	}

	std::function<double(double)> deviateFunction;
	std::function<double(double, const ParameterList&)> llFunction; //including the prior

	struct {
		double lower;
		double upper;
	} range;


	std::string type(void) const OVERRIDE {
		return "MH_Parameter";
	}

	AcceptanceTracker acceptanceTracker;
	AcceptanceTracker* getAcceptanceTracker(void) OVERRIDE {
		return &acceptanceTracker;
	}

protected:
	//friend class GibbsSampler;

	double _getNextSample(void) OVERRIDE;

};

class DecorrelatingStep : public GibbsParameter {
public:

	DecorrelatingStep(void) {
		_dimension = ParameterDimension::ZERO;
	}
	
	std::function<std::map<std::string, double>(const ParameterList&)> currentValuesFunction;
	std::function<std::map<std::string, double>(const ParameterList&)> candidateValuesFunction;
	std::function<double(const std::map<std::string, double>&, const ParameterList&)> llFunction;
	
	
	std::string type(void) const OVERRIDE {
		return "DecorrelatingStep";
	}

	AcceptanceTracker acceptanceTracker;
	AcceptanceTracker* getAcceptanceTracker(void) OVERRIDE {
		return &acceptanceTracker;
	}

protected:
	friend class GibbsSampler;

	double _getNextSample(void) OVERRIDE;
};


class VectorElement : public GibbsParameter {
public:

	VectorElement(void) {
		_dimension = ParameterDimension::SCALAR;
	}

	//neither of these is used, at present, but could be used for backtracking from this parameter to its creator
	//std::string sourceName;
	//unsigned int index;

	void update(void) OVERRIDE {
		//do nothing. this parameter is updated by the vector valued parameter that uses it
	}

	std::string type(void) const OVERRIDE {
		return "VectorElement";
	}

protected:

	double _getNextSample(void) OVERRIDE {
		return 0; //do nothing. Must override because it's pure virtual in the base class
	}

};


/*
This class does a block update of a vector-valued parameter.

In order to use this class, you must set all 4 properties:
startValues
ranges
deviateFunction
llFunction

Plus, before or after adding this parameter to the GibbsSampler, you must call createElements().
*/
class VectorMH_Parameter : public VectorParameter {
public:

	VectorMH_Parameter(void) {
		_dimension = ParameterDimension::VECTOR;
		startValues.resize(0); //whatev
	}

	VectorMH_Parameter(unsigned int size) {
		_dimension = ParameterDimension::VECTOR;
		startValues.resize(size, 0);
		ranges.resize(size, std::make_pair(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()));
	}

	std::vector<double> startValues;
	std::vector<std::pair<double, double>> ranges;

	std::function<std::vector<double>(const std::vector<double>&)> deviateFunction;
	std::function<double(std::vector<double>, const ParameterList&)> llFunction;

	

	void update(void) OVERRIDE;

	double value(void) const OVERRIDE {
		return std::numeric_limits<double>::quiet_NaN();
	}

	AcceptanceTracker acceptanceTracker;
	AcceptanceTracker* getAcceptanceTracker(void) OVERRIDE {
		return &acceptanceTracker;
	}

	std::string type(void) const OVERRIDE {
		return "VectorMH_Parameter";
	}

	//VectorParameterInterface functions
	void createElements(GibbsSampler* gibbs) OVERRIDE;

	const std::vector<std::string>& getElementNames(void) const OVERRIDE;

	std::vector<double> getCurrentValue(void) const OVERRIDE;

	std::vector<double> getIterationValue(unsigned int iteration) const OVERRIDE;

protected:

	double _getNextSample(void) OVERRIDE {
		return 0; // This is never called
	}

	void _parameterAddedToSampler(GibbsSampler* gs) OVERRIDE;

	void _parameterDeletedFromSampler(GibbsSampler* gs) OVERRIDE;


	std::vector<double> __getNextSamples(void);

};


class GibbsSampler {
public:

	GibbsSampler(void) :
		iterationsPerStatusUpdate(100),
		iterationCompleteCallback(nullptr),
		_storedIterations(0)
	{
		std::random_device rd;
		_generator.seed(rd());
	}

	void run(unsigned int samplesToCollect, bool clearExistingSamples = true);

	void clear(void);

	unsigned int iterationsPerStatusUpdate;

	//Takes a pointer to this and the current iteration. Returns true to continue, false to stop.
	std::function<bool(GibbsSampler*, unsigned int)> iterationCompleteCallback;

	void createConstantParameter(std::string name, std::string group, double value, std::string replaceName = "", bool removeAndAdd = true);
	void setParameterGroupToConstantValue(std::string group, double value);

	unsigned int getIndex(const std::string& parameterName) const;
	bool hasParameter(const std::string& parameterName) const;

	std::vector<GibbsParameter*> getParameters(void);
	std::vector<GibbsParameter*> getParameters(const std::vector<std::string>& param);

	std::vector<unsigned int> namesToIndices(std::vector<std::string> names) const;

	std::vector<double> getCurrentGroupValues(std::string group) const;


	void setCurrentParameterValue(std::string p, double v);
	const ParameterList& getCurrentParameterValues();
	ParameterList getParameterList(std::vector<std::string> parameterNames);
	ParameterList getParameterList(std::vector<unsigned int> parameterIndices);
	ParameterList getIterationParameterValues(unsigned int iteration);

	void updateDependentParameters(ParameterList* param) const;


	std::mt19937_64& getGenerator(void) {
		return _generator;
	}

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

		unsigned int index = _paramVector.size() - 1;

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
		std::vector<unsigned int>& gi = _groupToIndices[oldGroup];
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
	T* getParameter(unsigned int index) {
		GibbsParameter* p = _paramVector.at(index);
		return dynamic_cast<T*>(p);
	}

	template <typename T>
	const T* getParameter(unsigned int index) const {
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

	std::map<std::string, std::vector<unsigned int>> _groupToIndices;

	std::map<std::string, unsigned int> _nameToIndex;
	std::vector<GibbsParameter*> _paramVector;

	ParameterList _currentValues;

	unsigned int _storedIterations;

	std::mt19937_64 _generator;
	
	void _remakeIndices(void);

	void _makeDependentParameterList(void);
	std::vector<unsigned int> _continuouslyUpdatedDependentParameters;
	std::vector<unsigned int> _dependentParameters;
	std::vector<unsigned int> _independentParameters;

};

void outputDataVector(std::string filename, std::string header, std::vector<double> values);

