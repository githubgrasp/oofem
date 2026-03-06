#ifndef oofemcfg_h
#define oofemcfg_h

#include <cstddef>
#include "floatarrayf.h"

#ifdef _USE_SHARED
    #include "oofem_export.h"
#else
    #define OOFEM_EXPORT
    #define OOFEM_NO_EXPORT
#endif

OOFEM_EXPORT extern const char* PRG_VERSION;
OOFEM_EXPORT extern const char* PRG_VERSION_SHORT;
OOFEM_EXPORT extern const char* OOFEG_VERSION;
OOFEM_EXPORT extern const char* OOFEM_COPYRIGHT;
OOFEM_EXPORT extern const char* PRG_HEADER;
OOFEM_EXPORT extern const char* PRG_HEADER_SM;
OOFEM_EXPORT extern const char* HOST_TYPE;
OOFEM_EXPORT extern const char* HOST_NAME;
OOFEM_EXPORT extern const char* MODULE_LIST;
OOFEM_EXPORT extern const char* OOFEM_GIT_HASH;
OOFEM_EXPORT extern const char* OOFEM_GIT_REPOURL;
OOFEM_EXPORT extern const char* OOFEM_GIT_BRANCH;



namespace oofem {
template <std::size_t N> class FloatArrayF;
/** Define a type alias for 3D coordinates */
using Coordinates = FloatArrayF<3>;
} // end namespace oofem

#endif /* oofemcfg.h */
