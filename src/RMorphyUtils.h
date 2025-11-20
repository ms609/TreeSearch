#include <string.h>
#include "mpl.h"

#if defined (_WIN32) || defined(_WIN64) || defined(_WINDOWS)
#define strcasecmp _stricmp
#endif

MPLchtype R_mpl_str2chtype(const char *chtypename);
MPLgap_t R_mpl_str2gaptype(const char *chtypename);
