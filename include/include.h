#include "lapack.h"
#include "Matrix.h"
#include "Vector.h"

#include "RSPM.h"
#include "DIIS.h"

#include "Hub.h"

#define HDF5_STATUS_CHECK(status) if(status < 0) std::cerr << __FILE__ << ":" << __LINE__ << ": Problem with writing to file. Status code=" << status << std::endl;

/* vim: set ts=3 sw=3 expandtab :*/
