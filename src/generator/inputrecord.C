//
// file inputrecord.cc
//

#include "inputrecord.h"
#ifndef __MAKEDEPEND
#include <ctype.h>
#endif

InputRecord :: InputRecord()
{ }

InputRecord :: InputRecord(const InputRecord &src)
{ }



InputRecord &
InputRecord :: operator=(const InputRecord &src)
{
    return * this;
}

const char *
InputRecord :: strerror(IRResultType rt)
{
    switch ( rt ) {
    case IRRT_NOTFOUND:
        return "Missing Keyword"; // string literal is statically allocated, return safe

    case IRRT_BAD_FORMAT:
        return "Bad format";

    default:
        return "Unknown error";
    }
}

void
InputRecord :: report_error(const char *_class, const char *proc, const InputFieldType fieldID, const char *kwd,
                            IRResultType result, const char *file, int line)
{
  printf("error\n");
  exit(1);
  //   __OOFEM(file, line, "Input error: \"%s\", field keyword \"%s\" (fieldID=%d)\n%s::%s", strerror(result), kwd, fieldID, _class, proc);
}

