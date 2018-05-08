#ifndef SNUErrorHandler_H
#define SNUErrorHandler_H

// ROOT include(s):                                                                                                         
#include <TError.h>

/// Function printing log messages from ROOT                                                                                
extern void SNUErrorHandler( int level, Bool_t abort, const char* location,
                           const char* message );

#endif 



