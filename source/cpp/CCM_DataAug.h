
#include "GibbsSampler/GibbsSampler.h"

#include "CCM_Util.h"

namespace CatCont {

	struct daugConfig {

		unsigned int nCat;

		vector<string> indexedLabels; // TODO: Fill

	};

	// The parameters needed for a single observation
	struct daugParam_Obs {
		double label;
	};

	struct daugParam_Participant {

	};

	struct daugParam_Condition {

	};

	class daugModel {
	public:

		daugConfig dcfg;

		Data data;

		GibbsSampler gibbs;

		bool setup(daugCfg cfg, const Data& data);

	private:

		void _createParameters(void);

		double _singleObservationLL() const;

		// Functions with likelihood and prior
		double _obsLL_MH(double newLabel, const ParameterList& param, size_t partIndex, size_t condIndex, size_t obsIndex) const;

		

		double _labelDeviate(double labelIndex) const;

		double _labelNameToIndex(string name) const;
		string _labelIndexToName(double index) const;

	};




}

