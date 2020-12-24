/*---------------------------------------------------------------------------*\
  CFDEMcoupling - Open Source CFD-DEM coupling

  CFDEMcoupling is part of the CFDEMproject
  www.cfdem.com
                              Christoph Goniva, christoph.goniva@cfdem.com
                              Copyright 2009-2012 JKU Linz
                              Copyright 2012-     DCS Computing GmbH, Linz
------------------------------------------------------------------------------
License
  This file is part of CFDEMcoupling.

  CFDEMcoupling is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 3 of the License, or (at your
  option) any later version.

  CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with CFDEMcoupling; if not, write to the Free Software Foundation,
  Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
  cfdemCloud class managing DEM data for CFD-DEM coupling

Class
  Foam::cfdemCloud
\*---------------------------------------------------------------------------*/

#ifndef __CFDEM_CLOUD_H__
#define __CFDEM_CLOUD_H__

#include <memory>
#include <unordered_map>
#include <vector>

#include "IFstream.H"
#include "fvCFD.H"

#include "base/logging.h"
#include "base/memory/x_alloc.h"
#include "base/tensor/tensor.h"
#include "base/type_cast.h"
#include "cloud/cfdem_base.h"
#include "cloud/coupling_properties.h"
#include "cloud/of_version.h"
#include "cloud/particle_cloud.h"

// must behind include "of_version.h"
#if defined(version21) || defined(version16ext)
#include "turbulenceModel.H"
#elif defined(version15)
#include "RASModel.H"
#endif

namespace Foam {

// forward declarations
class liggghtsCommandModel;
class momCoupleModel;
class forceModel;
class averagingModel;
class voidFractionModel;
class dataExchangeModel;
class locateModel;

class cfdemCloud {
 public:
  //! \brief Constructor
  cfdemCloud(const fvMesh& mesh);

  //! \brief Destructor
  virtual ~cfdemCloud();

  //! \brief 重新分配内存
  virtual void reallocate();

  /*!
   * \brief 更新函数
   * \note used for cfdemSolverPiso
   * \param voidFraction  <[in, out] 小颗粒空隙率场
   * \param Us            <[in, out] 局部平均小颗粒速度场
   * \param U             <[in] 流体速度场
   */
  void evolve(volScalarField& voidFraction, volVectorField& Us, volVectorField& U);

  tmp<volScalarField> ddtVoidFraction() const;

 protected:
  /*!
   * \brief check if simulation is fully periodic
   * \return true if simulation is fully periodic
   */
  bool checkSimulationFullyPeriodic();

  /* --------------------------------- interfaces ---------------------------------------- */

 public:
  inline const fvMesh& mesh() const { return mesh_; }

  inline ParticleCloud& pCloud() { return parCloud_; }

  inline const CouplingProperties& cProps() const { return cProps_; }

  inline const IOdictionary& couplingPropertiesDict() const { return couplingPropertiesDict_; }

  inline const IOdictionary& liggghtsCommandsDict() const { return liggghtsCommandsDict_; }

  inline const std::vector<std::shared_ptr<liggghtsCommandModel>>& liggghtsCommandModels() const;

  inline const std::vector<std::shared_ptr<forceModel>>& forceModels() const;

  inline const std::vector<std::shared_ptr<momCoupleModel>>& momCoupleModels() const;

  inline std::shared_ptr<forceModel> forceM(int index) const;

  inline const dataExchangeModel& dataExchangeM() const;

  inline const voidFractionModel& voidFractionM() const;

  inline const locateModel& locateM() const;

  inline dataExchangeModel& dataExchangeM();

  inline voidFractionModel& voidFractionM();

  inline locateModel& locateM();

#if defined(version24Dev)
  inline const turbulenceModel& turbulence() const;
#elif defined(version21) || defined(version16ext)
#ifdef compre
  inline const compressible::turbulenceModel& turbulence() const;
#else
  inline const incompressible::turbulenceModel& turbulence() const;
#endif
#elif defined(version15)
  inline const incompressible::RASModel& turbulence() const;
#endif

  inline bool writeTimePassed() const { return writeTimePassed_; }

  /* ------------------------- interface of CouplingProperties --------------------------- */

  inline const std::vector<std::string>& forceModelList() const { return cProps_.forceModelList(); }

  inline const std::vector<std::string>& momCoupleModelList() const { return cProps_.momCoupleModelList(); }

  inline const std::vector<std::string>& liggghtsCommandModelList() const { return cProps_.liggghtsCommandModelList(); }

  inline bool verbose() const { return cProps_.verbose(); }

  inline bool solveFlow() const { return cProps_.solveFlow(); }

  inline const std::string& modelType() const { return cProps_.modelType(); }

  inline const std::string& turbulenceModelType() const { return cProps_.turbulenceModelType(); }

  inline bool allowUseSubCFDTimeStep() const { return cProps_.allowUseSubCFDTimeStep(); }

  inline int couplingInterval() const { return cProps_.couplingInterval(); }

  inline bool checkPeriodicCells() const { return cProps_.checkPeriodicCells(); }

  inline const Foam::vector& periodicCheckRange() const { return cProps_.periodicCheckRange(); }

  inline const std::string& ddtVoidFractionType() const { return cProps_.ddtVoidFractionType(); }

  /* ------------------------- interface of particleCloud ------------------------------- */

  inline int numberOfParticles() const { return parCloud_.numberOfParticles_; }

  inline bool numberOfParticlesChanged() const { return parCloud_.numberOfParticlesChanged_; }

  inline const base::CITensor1& particleOverMeshNumber() const { return parCloud_.particleOverMeshNumber_; }

  inline const base::CITensor1& findCellIDs() const { return parCloud_.findCellIDs_; }

  inline const base::CDTensor1& dimensionRatios() const { return parCloud_.dimensionRatios_; }

  inline const base::CDTensor1& volumes() const { return parCloud_.volumes_; }

  inline const base::CDTensor2& fAcc() const { return parCloud_.fAcc_; }

  inline const base::CDTensor2& impForces() const { return parCloud_.impForces_; }

  inline const base::CDTensor2& expForces() const { return parCloud_.expForces_; }

  inline const base::CDTensor2& particleWeights() const { return parCloud_.particleWeights_; }

  inline const base::CDTensor2& particleVolumes() const { return parCloud_.particleVolumes_; }

  inline const std::vector<base::CITensor1>& cellIDs() const { return parCloud_.cellIDs_; }

  inline const std::vector<base::CDTensor1>& voidFractions() const { return parCloud_.voidFractions_; }

  inline const std::vector<base::CDTensor1>& volumeFractions() const { return parCloud_.volumeFractions_; }

  inline double** radiiPtr() const { return parCloud_.radiiPtr_; }

  inline double** cdsPtr() const { return parCloud_.cdsPtr_; }

  inline double** positionsPtr() const { return parCloud_.positionsPtr_; }

  inline double** velocitiesPtr() const { return parCloud_.velocitiesPtr_; }

  inline double** angularVelocitiesPtr() const { return parCloud_.angularVelocitiesPtr_; }

  inline double** DEMForcesPtr() const { return parCloud_.DEMForcesPtr_; }

  inline double** DEMTorquesPtr() const { return parCloud_.DEMTorquesPtr_; }

  inline double** fluidVelPtr() const { return parCloud_.fluidVelPtr_; }

  inline const base::CDExTensor1& cds() const { return parCloud_.cds_; }

  inline const base::CDExTensor1& radii() const { return parCloud_.radii_; }

  inline const base::CDExTensor2& positions() const { return parCloud_.positions_; }

  inline const base::CDExTensor2& velocities() const { return parCloud_.velocities_; }

  inline const base::CDExTensor2& angularVelocities() const { return parCloud_.angularVelocities_; }

  inline const base::CDExTensor2& DEMForces() const { return parCloud_.DEMForces_; }

  inline const base::CDExTensor2& DEMTorques() const { return parCloud_.DEMTorques_; }

  inline const base::CDExTensor2& fluidVel() const { return parCloud_.fluidVel_; }

  inline double getRadius(int index) const { return parCloud_.getRadius(index); }

  inline Foam::vector getPosition(int index) const { return parCloud_.getPosition(index); }

  inline Foam::vector getVelocity(int index) const { return parCloud_.getVelocity(index); }

  inline Foam::vector getAngularVelocity(int index) const { return parCloud_.getAngularVelocity(index); }

  inline void setNumberOfParticles(int number) { parCloud_.setNumberOfParticles(number); }

 protected:
  //! \note 在当前类中一定要最先声明 couplingPropertiesDict_ 和 liggghtsCommandsDict_

  const fvMesh& mesh_;

  IOdictionary couplingPropertiesDict_;

  IOdictionary liggghtsCommandsDict_;

  CouplingProperties cProps_;

  ParticleCloud parCloud_;

  std::vector<std::shared_ptr<liggghtsCommandModel>> liggghtsCommandModels_;

  std::vector<std::shared_ptr<forceModel>> forceModels_;

  std::vector<std::shared_ptr<momCoupleModel>> momCoupleModels_;

  autoPtr<dataExchangeModel> dataExchangeModel_;

  autoPtr<voidFractionModel> voidFractionModel_;

  autoPtr<locateModel> locateModel_;

// autoPtr<averagingModel> averagingModel_;

// autoPtr<IOModel> IOModel_;

// autoPtr<probeModel> probeModel_;

// autoPtr<registryModel> registryModel_;

// autoPtr<clockModel> clockModel_;

// autoPtr<smoothingModel> smoothingModel_;

// autoPtr<meshMotionModel> meshMotionModel_;

#if defined(version24Dev)
  const turbulenceModel& turbulence_;
#elif defined(version21) || defined(version16ext)
#ifdef compre
  const compressible::turbulenceModel& turbulence_;
#else
  const incompressible::turbulenceModel& turbulence_;
#endif
#elif defined(version15)
  const incompressible::RASModel& turbulence_;
#endif

  //! \brief Multiphase Turbulence (e.g., slip-induced turbulence)
  volScalarField turbulenceMultiphase_;

  bool writeTimePassed_;

  bool arraysReallocated_;

  //! \brief used for ddt(voidFraction)
  volScalarField ddtVoidFraction_;

  /*!
   * \brief 是否隐式计算颗粒所受到的阻力
   *        true - 在每个耦合时间步，流体的速度和阻力系数都被传递到 DEM 中，从而在每个 DEM 时间步中，
   *               使用上一个耦合时间步中的阻力系数和流体速度，与当前颗粒速度一起计算颗粒受到的阻力
   *        false - 在每个耦合时间步中，流体对颗粒的阻力被传递到 DEM 中，并且在接下来的 DEM 时间步中，
   *                这个力保持不变，直到下个耦合时间步
   * \note default false
   */
  bool impDEMdrag_;

  /*!
   * \brief 在每个 DEM 时间步中，都将颗粒受到的力累计起来，然后在耦合时间步中传递给 CFD 计算，
   *        只有当 impDEMdrag_ = true 时才为真
   * \note default false
   */
  bool impDEMdragAcc_;

  //! \brief 显式、隐式分裂系数
  scalar impExpSplitFactor_;

  //! \brief Variable used to de-activate mirroring across periodic boundary conditions.
  Switch checkPeriodicCells_;

  /*!
   * \brief De-activation and tolerance variables, if set to (for a specific direction),
   *        the periodic check will NOT be done. Important for probing ambient points.
   *        Only read-in in case checkPeriodicCells is active.
   * \note default = (1, 1, 1), i.e., periodic checks will be done
   */
  vector wallPeriodicityCheckRange_;

  scalar wallPeriodicityCheckTolerance_;
};

}  // namespace Foam

#include "cloud/cfdem_cloud-inl.h"

#endif  // __CFDEM_CLOUD_H__
