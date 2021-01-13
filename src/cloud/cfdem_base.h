#ifndef __CFDEM_BASE_H__
#define __CFDEM_BASE_H__

#ifndef CFDEM_MIX_CLOUD
#define CFDEM_MIX_CLOUD 0
#endif

#include <type_traits>
#include "base/tensor/tensor.h"
#include "mpi.h"

namespace base {

void MPI_Barrier(int times = 0.1);

int numProc();

int procId();

template <int nDim, typename DType, typename Device, typename Alloc>
void MPI_Isend(const base::Tensor<nDim, DType, Device, Alloc>& tensor, const int destProc, int tag,
               MPI_Request* request) {
  if (std::is_same<double, DType>::value) {
    ::MPI_Isend(tensor.ptr(), tensor.mSize(), MPI_DOUBLE, destProc, tag, MPI_COMM_WORLD, request);
  } else if (std::is_same<int32_t, DType>::value) {
    ::MPI_Isend(tensor.ptr(), tensor.mSize(), MPI_INT, destProc, tag, MPI_COMM_WORLD, request);
  }
}

template <int nDim, typename DType, typename Device, typename Alloc>
void MPI_Irecv(const base::Tensor<nDim, DType, Device, Alloc>& tensor, const int srcPrco, int tag,
               MPI_Request* request) {
  if (std::is_same<double, DType>::value) {
    ::MPI_Irecv(tensor.ptr(), tensor.mSize(), MPI_DOUBLE, srcPrco, tag, MPI_COMM_WORLD, request);
  } else if (std::is_same<int32_t, DType>::value) {
    ::MPI_Irecv(tensor.ptr(), tensor.mSize(), MPI_INT, srcPrco, tag, MPI_COMM_WORLD, request);
  }
}

}  // namespace base

#endif  // __CFDEM_BASE_H__
