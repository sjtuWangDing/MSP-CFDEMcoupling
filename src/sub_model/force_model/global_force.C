#include "./global_force.h"
#include "sub_model/void_fraction_model/void_fraction_model.h"

namespace Foam {

globalForce::globalForce(cfdemCloud& cloud)
    : cloud_(cloud),
      impParticleForce_(IOobject("impParticleForce", cloud.mesh().time().timeName(), cloud.mesh(),
                                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                        cloud.mesh(), dimensionedVector("zero", dimensionSet(1, 1, -2, 0, 0),
                                                        vector(0, 0, 0))  // [N] == [kg * m / s^2]
                        ),
      expParticleForce_(IOobject("expParticleForce", cloud.mesh().time().timeName(), cloud.mesh(),
                                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                        cloud.mesh(), dimensionedVector("zero", dimensionSet(1, 1, -2, 0, 0),
                                                        vector(0, 0, 0))  // [N] == [kg * m / s^2]
                        ) {}

globalForce::~globalForce() {}

void globalForce::updateExpandedCellMap() {
  // 清空数据
  expandedCellMap_.clear();
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkMiddleParticle(index)) {
      double radius = cloud_.getRadius(index);                       // 颗粒半径
      Foam::vector particlePos = cloud_.getPosition(index);          // 颗粒中心坐标
      int findExpandedCellID = cloud_.findExpandedCellIDs()[index];  // 扩展网格ID
      // 插入颗粒扩展网格的集合
      expandedCellMap_.insert(std::make_pair(index, std::unordered_set<int>()));
      if (findExpandedCellID >= 0) {
        // 构建扩展网格集合
        cloud_.voidFractionM().buildExpandedCellSet(expandedCellMap_[index], findExpandedCellID, particlePos, radius,
                                                    cloud_.expandedCellScale());
      }
    }
  }
  base::MPI_Barrier();
}

}  // namespace Foam
