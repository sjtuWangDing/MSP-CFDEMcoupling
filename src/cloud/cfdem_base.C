#include <unistd.h>
#include "cloud/cfdem_base.h"

namespace base {

void MPI_Barrier(int times) {
  usleep(times * 1000000);
  MPI_Barrier(MPI_COMM_WORLD);
}

int numProc() {
  int numProc = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  return numProc;
}

int procId() {
  int procId = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  return procId;
}

}  // namespace base
