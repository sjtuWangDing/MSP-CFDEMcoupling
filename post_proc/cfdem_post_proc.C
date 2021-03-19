#include "./cfdem_post_proc.h"

int main() {
  std::string dumFileName("dump_fineParticles_vm.run");
  cfdemPostProc cfdemPostProc;
  double DEMDt = 0.000002;

  std::vector<Node> NodeVec = cfdemPostProc.ReadDmpRunFile(dumFileName, DEMDt);

#if 0
  std::vector<std::pair<double, Foam::vector> > volVec1 =
      cfdemPostProc.volAverage(Foam::vector(0.0, 0.0, 0.041), Foam::vector(0.005, 0.00143, 0.0194), 12, NodeVec);
  Foam::vector sumVel = Foam::vector::zero;
  int number = 0;
  for (const auto& pair : volVec1) {
    std::cout << (pair.first - 11) << " " << pair.second[0] << " " << pair.second[1] << " " << pair.second[2]
              << std::endl;
    sumVel += pair.second;
    number += 1;
  }
  Foam::vector aveVel = sumVel / number;
  std::cout << aveVel[0] << " " << aveVel[1] << " " << aveVel[2] << " " << std::endl;
#elif 1
  double dimX = 0.0025, dimY = 0.001, dimZ = 0.003, dZ = 0.002;
  Foam::vector rectDim(dimX, dimY, dimZ);
  Foam::vector startPoint = Foam::vector::zero;

  startPoint = Foam::vector(0.0, 0.0, 0.0);
  int nrow = 1, ncol = 1;
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      startPoint[0] = i * 0.001;
      startPoint[1] = j * 0.001;
      startPoint[2] = 0.0;
      std::cout << "start point: " << startPoint[0] << " " << startPoint[1] << " " << startPoint[2] << std::endl;
      for (double height = 0.0; height < 0.106 + Foam::SMALL; height += dZ) {
        Foam::vector timeAveVel = cfdemPostProc.timeAverage(startPoint, rectDim, 12, DEMDt, NodeVec);
        std::cout << startPoint[2] << " " << timeAveVel[0] << " " << timeAveVel[1] << " " << timeAveVel[2] << std::endl;
        startPoint.z() = startPoint.z() + dZ;
      }
    }
  }
#else
  double radius = 0.001, dimZ = 0.004, dZ = 0.002;
  Foam::vector startPoint(0.001, 0.0, 0.0);
  for (double height = 0.0; height < 0.104 + Foam::SMALL; height += dZ) {
    Foam::vector timeAveVel = cfdemPostProc.timeAverage(startPoint, radius, dimZ, 11, DEMDt, NodeVec);
    std::cout << startPoint[2] << " " << timeAveVel[0] << " " << timeAveVel[1] << " " << timeAveVel[2] << std::endl;
    startPoint.z() = startPoint.z() + dZ;
  }
#endif
  return 0;
}
