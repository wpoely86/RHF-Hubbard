#include "lapack.h"
#include "Matrix.h"
#include "Vector.h"

#include "LibInt.h"

#include "R.h"
#include "Gauss.h"
#include "input.h"

#include "preamble.h"

#include "CI_SPM.h"
#include "CI_TPM.h"

#include "CI_SPPM.h"
#include "CI_SPPM_m.h"

#include "CI_TPPM.h"

#include "CartInt.h"

#include "Transform.h"

#include "SI_SPM.h"
#include "SI_TPM.h"

#include "SphInt.h"

#include "Tools.h"

#include "RSPM.h"
#include "DIIS.h"

#define HDF5_STATUS_CHECK(status) if(status < 0) std::cerr << __FILE__ << ":" << __LINE__ << ": Problem with writing to file. Status code=" << status << std::endl;

/* vim: set ts=3 sw=3 expandtab :*/
