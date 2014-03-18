#ifndef RAPIDJSON_COOLPROP_H
#define RAPIDJSON_COOLPROP_H

// On PowerPC, we are going to use the stdint.h integer types and not let rapidjson use its own
#if defined(__powerpc__)
typedef unsigned int UINT32;
#include "stdint.h"
#define RAPIDJSON_NO_INT64DEFINE
#endif

#include "rapidjson/rapidjson.h"
#include "rapidjson/document.h"
#include "rapidjson/filestream.h"	// wrapper of C stream for prettywriter as output
#include "rapidjson/prettywriter.h"	// for stringify JSON
#include "rapidjson/stringbuffer.h" // for string buffer

#endif
