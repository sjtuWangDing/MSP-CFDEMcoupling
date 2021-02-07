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
#include "base/run_time_selection_tables.h"
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

namespace base {

//! \brief multi parameter
template <typename... Types>
inline void MPI_Info(const Types&... args, bool verbose) {
  if (verbose) {
    oStreamPrint(Foam::Info, args...);
    Foam::Info << Foam::endl;
  }
  MPI_Barrier();
}

//! \brief single parameter
template <typename Type>
inline void MPI_Info(const Type& arg, bool verbose) {
  if (verbose) {
    oStreamPrint(Foam::Info, arg);
    Foam::Info << Foam::endl;
  }
  MPI_Barrier();
}

//! \brief multi parameter for const char *
template <typename... CTypes>
inline void MPI_Info(const CTypes*... args, bool verbose) {
  if (verbose) {
    oStreamPrint(Foam::Info, args...);
    Foam::Info << Foam::endl;
  }
  MPI_Barrier();
}

//! \brief single parameter for const char *
template <typename CType>
inline void MPI_Info(const CType* arg, bool verbose) {
  if (verbose) {
    oStreamPrint(Foam::Info, arg);
    Foam::Info << Foam::endl;
  }
  MPI_Barrier();
}

}  // namespace base

namespace Foam {

// forward declarations
class liggghtsCommandModel;
class momCoupleModel;
class forceModel;
class globalForce;
class voidFractionModel;
class dataExchangeModel;
class locateModel;
class averagingModel;

typedef std::unordered_map<std::string, std::shared_ptr<momCoupleModel>> momCoupleModelMap;

class cfdemCloud {
 public:
  //! \brief Constructor
  cfdemCloud(const fvMesh& mesh);

  //! \brief Destructor
  virtual ~cfdemCloud();

  //! \brief Runtime type information
  cfdemTypeName("cfdemCloud");

  virtual double expandedCellScale() const {
    FatalError << __func__
               << " not implement in Foam::cfdemCloud, please use Foam::cfdemCloudSemi or Foam::cfdemCloudMix"
               << abort(FatalError);
    return 0.0;
  }

  virtual bool checkFineParticle(int index) const { return true; }

  virtual bool checkMiddleParticle(int index) const { return true; }

  virtual bool checkCoarseParticle(int index) const { return true; }

  /*!
   * \brief 更新函数
   * \note used for cfdemSolverPiso
   * \param U      <[in] 流体速度场
   * \param voidF  <[in, out] 小颗粒空隙率场
   * \param Us     <[in, out] 局部平均小颗粒速度场
   * \param Ksl    <[in, out] 动量交换场
   */
  virtual void evolve(volVectorField& U, volScalarField& voidF, volVectorField& Us, volScalarField& Ksl);

  void calcDDtVoidFraction(volScalarField& voidFraction);

  tmp<volScalarField> ddtVoidFraction() const;

  tmp<volScalarField> voidFractionNuEff(volScalarField& voidFraction) const;

  tmp<fvVectorMatrix> divVoidFractionTau(volVectorField& U, volScalarField& voidFraction) const;

  template <typename Field>
  inline void scaleWithVcell(Field& field) const;

 protected:
  //! \brief 初始化函数，用于在构造函数中执行
  void init();

  //! \brief 重新分配内存
  void reallocate();

  void getDEMData();

  void giveDEMData() const;

  void printParticleInfo() const;

  //! \brief reset field
  void resetField();

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

  inline const std::vector<std::shared_ptr<liggghtsCommandModel>>& liggghtsCommandModels() const {
    return liggghtsCommandModels_;
  }

  inline const std::vector<std::shared_ptr<forceModel>>& forceModels() const { return forceModels_; }

  inline const momCoupleModelMap& momCoupleModels() const { return momCoupleModels_; }

  const globalForce& globalF() const;

  const dataExchangeModel& dataExchangeM() const;

  const voidFractionModel& voidFractionM() const;

  const locateModel& locateM() const;

  const averagingModel& averagingM() const;

  globalForce& globalF();

  dataExchangeModel& dataExchangeM();

  voidFractionModel& voidFractionM();

  locateModel& locateM();

  averagingModel& averagingM();

#if defined(version24Dev)
  const turbulenceModel& turbulence() const;
#elif defined(version21) || defined(version16ext)
#ifdef compre
  const compressible::turbulenceModel& turbulence() const;
#else
  const incompressible::turbulenceModel& turbulence() const;
#endif
#elif defined(version15)
  inline const incompressible::RASModel& turbulence() const;
#endif

  inline bool writeTimePassed() const { return writeTimePassed_; }

  inline bool meshHasUpdated() const { return meshHasUpdated_; }

  inline void setMeshHasUpdated(bool meshHasUpdated) { meshHasUpdated_ = meshHasUpdated; }

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

  inline bool useGuoBBOEquation() const { return cProps_.useGuoBBOEquation(); }

  inline bool useDDtVoidFraction() const { return cProps_.useDDtVoidFraction(); }

  inline double refineMeshSkin() const { return cProps_.refineMeshSkin(); }

  inline double minCoarseParticleRadius() const { return cProps_.minCoarseParticleRadius(); }

  /* ------------------------- interface of particleCloud ------------------------------- */

  inline int numberOfParticles() const { return parCloud_.numberOfParticles_; }

  inline bool numberOfParticlesChanged() const { return parCloud_.numberOfParticlesChanged_; }

  inline const base::CITensor1& particleOverMeshNumber() const { return parCloud_.particleOverMeshNumber_; }

  inline const base::CITensor1& particleRootProcIDs() { return parCloud_.particleRootProcIDs_; }

  inline const base::CITensor1& findCellIDs() { return parCloud_.findCellIDs_; }

  inline const base::CITensor1& findMpiCellIDs() { return parCloud_.findMpiCellIDs_; }

  inline const base::CITensor1& findExpandedCellIDs() { return parCloud_.findExpandedCellIDs_; }

  inline const base::CDTensor1& dimensionRatios() const { return parCloud_.dimensionRatios_; }

  inline const base::CDTensor2& impForces() const { return parCloud_.impForces_; }

  inline const base::CDTensor2& expForces() const { return parCloud_.expForces_; }

  inline const std::vector<base::CITensor1>& cellIDs() const { return parCloud_.cellIDs_; }

  inline const std::vector<base::CDTensor1>& voidFractions() const { return parCloud_.voidFractions_; }

  inline const std::vector<base::CDTensor1>& volumeFractions() const { return parCloud_.volumeFractions_; }

  inline const std::vector<base::CDTensor1>& particleWeights() const { return parCloud_.particleWeights_; }

  inline const std::vector<base::CDTensor1>& particleVolumes() const { return parCloud_.particleVolumes_; }

  inline double** radiiPtr() const { return parCloud_.radiiPtr_; }

  inline double** cdsPtr() const { return parCloud_.cdsPtr_; }

  inline double** positionsPtr() const { return parCloud_.positionsPtr_; }

  inline double** velocitiesPtr() const { return parCloud_.velocitiesPtr_; }

  inline double** initVelocitiesPtr() const { return parCloud_.initVelocitiesPtr_; }

  inline double** angularVelocitiesPtr() const { return parCloud_.angularVelocitiesPtr_; }

  inline double** DEMForcesPtr() const { return parCloud_.DEMForcesPtr_; }

  inline double** DEMTorquesPtr() const { return parCloud_.DEMTorquesPtr_; }

  inline double** fluidVelPtr() const { return parCloud_.fluidVelPtr_; }

  inline const base::CDExTensor1& cds() const { return parCloud_.cds_; }

  inline const base::CDExTensor1& radii() const { return parCloud_.radii_; }

  inline const base::CDExTensor2& positions() const { return parCloud_.positions_; }

  inline const base::CDExTensor2& velocities() const { return parCloud_.velocities_; }

  inline const base::CDExTensor2& initVelocities() const { return parCloud_.initVelocities_; }

  inline const base::CDExTensor2& angularVelocities() const { return parCloud_.angularVelocities_; }

  inline const base::CDExTensor2& DEMForces() const { return parCloud_.DEMForces_; }

  inline const base::CDExTensor2& DEMTorques() const { return parCloud_.DEMTorques_; }

  inline const base::CDExTensor2& fluidVel() const { return parCloud_.fluidVel_; }

  inline double getRadius(int index) const { return parCloud_.getRadius(index); }

  inline double getDimensionRatio(int index) const { return parCloud_.getDimensionRatio(index); }

  inline Foam::vector getPosition(int index) const { return parCloud_.getPosition(index); }

  inline Foam::vector getVelocity(int index) const { return parCloud_.getVelocity(index); }

  inline Foam::vector getInitVelocity(int index) const { return parCloud_.getInitVelocity(index); }

  inline Foam::vector getAngularVelocity(int index) const { return parCloud_.getAngularVelocity(index); }

  inline void setNumberOfParticles(int number) {
    Info << "number of particles: " << number << endl;
    parCloud_.setNumberOfParticles(number);
  }

 protected:
  //! \note 在当前类中一定要最先声明 couplingPropertiesDict_ 和 liggghtsCommandsDict_

  const fvMesh& mesh_;

  IOdictionary couplingPropertiesDict_;

  IOdictionary liggghtsCommandsDict_;

  CouplingProperties cProps_;

  ParticleCloud parCloud_;

  bool writeTimePassed_;

  /*!
   * \brief 判断 mesh 是否被更新过
   * \note 在求解器中使用 dynamic mesh，如果 mesh 更新，则 mesh.update() 返回 true，否则返回 false
   */
  bool meshHasUpdated_;

  bool validCouplingStep_;

  std::vector<std::shared_ptr<liggghtsCommandModel>> liggghtsCommandModels_;

  std::vector<std::shared_ptr<forceModel>> forceModels_;

  momCoupleModelMap momCoupleModels_;

  autoPtr<globalForce> globalForce_;

  autoPtr<dataExchangeModel> dataExchangeModel_;

  autoPtr<voidFractionModel> voidFractionModel_;

  autoPtr<locateModel> locateModel_;

  autoPtr<averagingModel> averagingModel_;

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

  //! \brief used for ddt(voidFraction)
  volScalarField ddtVoidFraction_;

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

template <typename Field>
void cfdemCloud::scaleWithVcell(Field& field) const {
  forAll(field, cellID) { field[cellID] *= mesh_.V()[cellID]; }
}

}  // namespace Foam

#endif  // __CFDEM_CLOUD_H__
