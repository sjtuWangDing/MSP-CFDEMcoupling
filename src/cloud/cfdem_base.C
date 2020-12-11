#include <unistd.h>
#include "mpi.h"

namespace base {

void MPI_Barrier(int times) {
  usleep(times * 1000000);
  MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace base
