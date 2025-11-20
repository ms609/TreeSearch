#include <Rcpp.h>
using namespace Rcpp;

#include "RMorphyUtils.h"

MPLchtype R_mpl_str2chtype(const char *chtypename)
{

	if(!strcasecmp(chtypename, "fitch")){
		return FITCH_T;
	}
	else if(!strcasecmp(chtypename, "wagner")){
		return WAGNER_T;
	}
	else if(!strcasecmp(chtypename, "dollo")){
		return DOLLO_T;
	}
	else if(!strcasecmp(chtypename, "irreversible")){
		return IRREVERSIBLE_T;
	}
	else if(!strcasecmp(chtypename, "user")){
		return USERTYPE_T;
	}

	return MAX_CTYPE;
}

MPLgap_t R_mpl_str2gaptype(const char *chtypename)
{

	if(!strcasecmp(chtypename, "inapplicable")){
		return GAP_INAPPLIC;
	}
	else if(!strcasecmp(chtypename, "missing")){
		return GAP_MISSING;
	}
	else if(!strcasecmp(chtypename, "newstate")){
		return GAP_NEWSTATE;
	}

	return GAP_MAX;
}
