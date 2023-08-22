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
#include "UtilityFunctions.h"

class GibbsSampler; //Forward declaration

typedef std::map<std::string, double> ParameterList; //A map string -> double is a conventient way to pass around parameter values

// Simple way to pass named values
typedef std::map<std::string, double> ParamMap;

// Idea: Interface?
struct ParamMapInterface {
	virtual double get(const std::string& name) const = 0;
};

// Complex way to pass named values, with optional indexing for speed.
class ParamContainer {
public:

	static const size_t BadIndex = -1;

	struct ParamProxy {

		bool setup(const ParamContainer& source, const std::string& name) {
			this->name = name;
			this->index = source.getIndex(name);
			return this->index != ParamContainer::BadIndex;
		}

		std::string name;
		size_t index = ParamContainer::BadIndex;
		double value = 0;

		double getValue(const ParamContainer& source) const {
			return source.get(this->index);
		}

		void updateValue(const ParamContainer& source) {
			this->value = this->getValue(source);
		}
		
	};

	ParamProxy makeProxy(const std::string& name) const {
		ParamProxy rval;
		rval.setup(*this, name);
		return rval;
	}

	bool setup(const std::vector<std::string>& names) {
		std::vector<double> values(names.size(), 0);
		return setup(names, values);
	}

	bool setup(const std::vector<std::string>& names, const std::vector<double>& values) {
		if (names.size() != values.size()) {
			return false;
		}

		//_names = names;
		_values = values;
		for (size_t i = 0; i < names.size(); i++) {
			_nameToIndex[names[i]] = i;
		}

		return true;
	}

	bool setValues(const std::vector<double>& values) {
		if (values.size() == this->_values.size()) {
			this->_values = values;
			return true;
		}
		return false;
	}

	double& operator[](const std::string& name) {
		// Check if need to create parameter
		if (!this->hasName(name)) {
			_nameToIndex[name] = _values.size();
			_values.push_back(0);
		}

		size_t ind = _nameToIndex.at(name);
		return _values[ind];
	}

	const double& operator[](const std::string& name) const {
		return this->at(name);
	}

	const double& at(const std::string& name) const {
		// No creating parameters in this function

		size_t ind = _nameToIndex.at(name);
		return _values[ind];
	}

	void set(const std::string& name, double value) {
		size_t ind = _nameToIndex.at(name);
		this->set(ind, value);
	}

	double get(const std::string& name) const {
		size_t ind = _nameToIndex.at(name);
		return this->get(ind);
	}

	void set(size_t ind, double value) {
		_values[ind] = value;
	}

	double get(size_t ind) const {
		return _values[ind];
	}

	void set(ParamProxy pp, double value) {
		_values[pp.index] = value;
	}

	double get(ParamProxy pp) const {
		return _values[pp.index];
	}

	bool hasName(const std::string& name) const {
		return _nameToIndex.find(name) != _nameToIndex.end();
	}

	size_t getIndex(const std::string& name) const {
		const auto& it = _nameToIndex.find(name);
		if (it == _nameToIndex.end()) {
			return BadIndex;
		}
		return it->second;
	}

	std::vector<std::string> reconstructOrderedNames(void) const {
		std::vector<std::string> rval(_values.size());
		for (const auto& it : _nameToIndex) {
			rval[it.second] = it.first;
		}
		return rval;
	}

private:

	//std::vector<std::string> _names;
	std::vector<double> _values;
	std::map<std::string, size_t> _nameToIndex; // Index in this container

};


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

	// This function is unsafe to use
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


class SectionMarker : public GibbsParameter {

	SectionMarker(void) :
		GibbsParameter()
	{
		group = this->type();
		_dimension = ParameterDimension::ZERO;
	}

	virtual std::string type(void) const override {
		return "SectionMarker";
	}

	std::function<void(std::string)> updateCallback; // argument is name

	virtual void update(void) override {
		// do stuff
		updateCallback(this->name);
	}

};


class BaseVectorParameter : public GibbsParameter {
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

	//std::function<double(const ParameterList&)> samplingFunction;
	std::function<double(const ParamContainer&)> samplingFunction;

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

	bool updateContinuously; // If false, this parameter value is only updated once per iteration, not when other parameters are updated.

	bool hasDependency(void) OVERRIDE {
		return true;
	}

	// Names of parameters on which this parameter depends. 
	// Returning an empty vector is treated as dependent on all other parameters (i.e. no parameter may depend on this parameter).
	std::vector<std::string> getDependencyNames(void) const {
		return {};
	}

	//virtual double evaluate(const ParameterList& param) const = 0;
	virtual double evaluate(const ParamContainer& param) const = 0;
};

class DependentParameter : public GenericDependentParameter {
public:

	DependentParameter(void) {
		_dimension = ParameterDimension::SCALAR;
	}

	DependentParameter(std::string name_, std::string group_, std::string source) {
		_dimension = ParameterDimension::SCALAR;

		this->name = name_;
		this->group = group_;

		sourceParameter = source;
	}

	std::string sourceParameter; //The name of the parameter from which to take its value

	std::string type(void) const OVERRIDE {
		return "DependentParameter";
	}

	//double evaluate(const ParameterList& param) const OVERRIDE;
	double evaluate(const ParamContainer& param) const OVERRIDE;


	double value(void) const OVERRIDE;

	double getSample(unsigned int iteration) const OVERRIDE;
	std::vector<double>& getSamples(void) OVERRIDE;

	virtual bool hasDependency(void) OVERRIDE {
		return true;
	}

	virtual std::vector<std::string> getDependencyNames(void) const OVERRIDE {
		return { this->sourceParameter };
	}

protected:

	double _getNextSample(void) OVERRIDE;

};

// This class has a problem with specifying which parameters it depends on.
class CalculatedParameter : public GenericDependentParameter {
public:
	CalculatedParameter(void) {
		_dimension = ParameterDimension::SCALAR;
	}

	//void setup(std::function<double(const ParamContainer&)> samplingFunction_, const std::vector<std::string>& sourceNames_);

	//An arbitrary function that can do whatever (but probably calculate the parameter based on other parameters).
	//std::function<double(const ParameterList&)> samplingFunction;
	std::function<double(const ParamContainer&)> samplingFunction;

	std::vector<std::string> sourceParameters; // Optional but helpful

	std::string type(void) const OVERRIDE {
		return "CalculatedParameter";
	}

	//double evaluate(const ParameterList& param) const OVERRIDE;
	double evaluate(const ParamContainer& param) const OVERRIDE;

	double value(void) const OVERRIDE;

	double getSample(unsigned int iteration) const OVERRIDE;
	std::vector<double>& getSamples(void) OVERRIDE;

	virtual bool hasDependency(void) OVERRIDE {
		return true;
	}

	virtual std::vector<std::string> getDependencyNames(void) const OVERRIDE {
		return this->sourceParameters;
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
class ConstantVectorParameter : public BaseVectorParameter {
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

/* Metropolis-Hastings step parameter.

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
	//std::function<double(double, const ParameterList&)> llFunction; //including the prior
	std::function<double(double, const ParamContainer&)> llFunction; //including the prior

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

	//std::function<std::map<std::string, double>(const ParameterList&)> currentValuesFunction;
	//std::function<std::map<std::string, double>(const ParameterList&)> candidateValuesFunction;
	//std::function<double(const std::map<std::string, double>&, const ParameterList&)> llFunction;

	std::function<ParamMap(const ParamContainer&)> currentValuesFunction;
	std::function<ParamMap(const ParamContainer&)> candidateValuesFunction;
	std::function<double(const ParamMap&, const ParamContainer&)> llFunction;


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
class VectorMH_Parameter : public BaseVectorParameter {
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

	//std::function<double(std::vector<double>, const ParameterList&)> llFunction;
	std::function<double(std::vector<double>, const ParamContainer&)> llFunction;



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



