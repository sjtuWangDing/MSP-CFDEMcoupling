#ifndef __CFDEM_BASE_H__
#define __CFDEM_BASE_H__

#ifndef CFDEM_USE_TENSOR
#define CFDEM_USE_TENSOR 1
#endif

#ifndef CFDEM_MIX_CLOUD
#define CFDEM_MIX_CLOUD 0
#endif

namespace base {

void MPI_Barrier(int times = 0.1);

}  // namespace base

#endif  // __CFDEM_BASE_H__
