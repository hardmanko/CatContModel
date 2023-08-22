#include "CCM_Database.h"

namespace CatCont {


	void CCM_Database::setPartParam(PartInd part, ParamInd param, double val) {

		ParticipantParameters& pp = _parts[part].partParam;

		switch (param) {
		case ParamInd::pMem:
			pp.pMem = val;
		}

		//_parts[part].partParam.

	}

} // namespace CatCont