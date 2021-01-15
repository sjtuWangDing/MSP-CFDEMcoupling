#include "./global_force.h"
#include "sub_model/void_fraction_model/void_fraction_model.h"

namespace Foam {

globalForce::globalForce(cfdemCloud& cloud)
    : cloud_(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict("globalForceProps")),
      gravityFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("gravityFieldName", "g").c_str()),
      densityFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("densityFieldName", "rho").c_str()),
      velFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("velFieldName", "U").c_str()),
      voidFractionFieldName_(
          subPropsDict_.lookupOrDefault<Foam::word>("voidFractionFieldName", "voidFraction").c_str()),
      phiFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("phiFieldName", "phi").c_str()),
#if defined(version21)
      g_(cloud.mesh().lookupObject<uniformDimensionedVectorField>(gravityFieldName_)),
#elif defined(version16ext) || defined(version15)
      g_(dimensionedVector(
             cloud.mesh().lookupObject<IOdictionary>("environmentalProperties").lookup(environmentalProperties))
             .value()),
#endif
      rho_(cloud.mesh().lookupObject<volScalarField>(densityFieldName_)),
      U_(cloud.mesh().lookupObject<volVectorField>(velFieldName_)),
      voidFraction_(cloud.mesh().lookupObject<volScalarField>(voidFractionFieldName_)),
      phi_(cloud.mesh().lookupObject<surfaceScalarField>(phiFieldName_)),
      ddtU_(IOobject("ddtU", cloud.mesh().time().timeName(), cloud.mesh(), IOobject::READ_IF_PRESENT,
                     IOobject::AUTO_WRITE),
            cloud.mesh(),
            dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0), vector(0, 0, 0))),  // [ddtU] == [m / s^2])
      impParticleForce_(IOobject("impParticleForce", cloud.mesh().time().timeName(), cloud.mesh(),
                                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                        cloud.mesh(), dimensionedVector("zero", dimensionSet(1, 1, -2, 0, 0),
                                                        vector(0, 0, 0))),  // [N] == [kg * m / s^2]
      expParticleForce_(IOobject("expParticleForce", cloud.mesh().time().timeName(), cloud.mesh(),
                                 IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                        cloud.mesh(), dimensionedVector("zero", dimensionSet(1, 1, -2, 0, 0),
                                                        vector(0, 0, 0)))  // [N] == [kg * m / s^2]
{
}

globalForce::~globalForce() {}

//! \brief 更新颗粒覆盖的扩展网格索引
void globalForce::updateKernel() {
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

//! \brief 每一次耦合中，在 set force 前执行
void globalForce::initBeforeSetForce() {
  // init data
  expandedCellMap_.clear();
  backgroundUfluidMap_.clear();
  backgroundVoidFractionMap_.clear();
  backgroundDDtUMap_.clear();
  // 计算 ddtU field
  ddtU_ = fvc::ddt(U_) + fvc::div(phi_, U_);
  // reset data
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkMiddleParticle(index)) {
      double radius = cloud_.getRadius(index);                       // 颗粒半径
      Foam::vector particlePos = cloud_.getPosition(index);          // 颗粒中心坐标
      int findExpandedCellID = cloud_.findExpandedCellIDs()[index];  // 扩展网格ID
      // 计算颗粒覆盖的扩展网格集合
      expandedCellMap_.insert(std::make_pair(index, std::unordered_set<int>()));
      if (findExpandedCellID >= 0) {
        cloud_.voidFractionM().buildExpandedCellSet(expandedCellMap_[index], findExpandedCellID, particlePos, radius,
                                                    cloud_.expandedCellScale());
      }
      // 计算背景流体速度
      backgroundUfluidMap_.insert(std::make_pair(index, getBackgroundFieldValue(index, U_)));
      // 计算背景流体空隙率
      backgroundVoidFractionMap_.insert(
          std::make_pair(index, getBackgroundFieldValue<false, 1, volScalarField, scalar>(index, voidFraction_)));
      // 计算背景流体的 ddtU
      backgroundDDtUMap_.insert(std::make_pair(index, getBackgroundFieldValue(index, ddtU_)));
      // 计算初始时刻的相对速度(no need to clear initUrMap_ before set)
      if (cloud_.dataExchangeM().isFirstCouplingStep()) {
        initUrMap_.insert(std::make_pair(index, backgroundUfluidMap_[index] - cloud_.getVelocity(index)));
      }
    }
  }
}

//! \brief 每一次耦合中，在 set force 后执行
void globalForce::endAfterSetForce() {
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // (1) 更新颗粒速度，所有处理器都必须更新，当颗粒从一个计算域运动到另一个计算域，新的计算域必须有颗粒的 prevU
    if (prevParticleVelMap_.end() == prevParticleVelMap_.find(index)) {
      prevParticleVelMap_.insert(std::make_pair(index, cloud_.getVelocity(index)));
    } else {
      prevParticleVelMap_[index] = cloud_.getVelocity(index);
    }
    // (2) 更新 ddtUrHistoryMap_，所有处理器都必须更新
    if (cloud_.checkMiddleParticle(index)) {
      int procId = base::procId();                         // 处理器编号
      int numProc = base::numProc();                       // 处理器数量
      int rootProc = cloud_.particleRootProcIDs()[index];  // 主节点编号，即颗粒所在的处理器编号
      // 只有一个节点直接跳过
      if (1 == numProc) {
        continue;
      }
      std::vector<Foam::vector>& ddtUrHistoryVec = getDDtUrHistory(index);
      // data buffer
      std::vector<double> ddtUr(3, 0.0);
      // 颗粒中心所在的主节点作为广播节点填充 data buffer
      if (rootProc == procId) {
        CHECK(!ddtUrHistoryVec.empty()) << __func__ << ": ddtUrHistoryVec of particle " << index
                                        << " in the root proc is empty";
        int size = ddtUrHistoryVec.size();
        ddtUr[0] = ddtUrHistoryVec[size - 1][0];
        ddtUr[1] = ddtUrHistoryVec[size - 1][1];
        ddtUr[2] = ddtUrHistoryVec[size - 1][2];
      }
      MPI_Bcast(ddtUr.data(), ddtUr.size(), MPI_DOUBLE, rootProc, MPI_COMM_WORLD);
      // 所有的非主节点保存数据
      if (rootProc != procId) {
        ddtUrHistoryVec.push_back(Foam::vector(ddtUr[0], ddtUr[1], ddtUr[2]));
      }
    }
    base::MPI_Barrier();
  }
}

}  // namespace Foam
