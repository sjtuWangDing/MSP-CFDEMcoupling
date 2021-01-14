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
  This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
  and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

Class
  Foam::CouplingProperties
\*---------------------------------------------------------------------------*/

#ifndef __COUPLING_PROPERTIES_H__
#define __COUPLING_PROPERTIES_H__

#include <vector>
#include "cloud/cfdem_base.h"
#include "dictionary.H"
#include "fvMesh.H"
#include "vector.H"

namespace Foam {

class CouplingProperties {
 public:
  //! \brief Constructor
  CouplingProperties(const fvMesh& mesh, const IOdictionary& couplingPropertiesDict,
                     const IOdictionary& liggghtsCommandsDict);

  /* --------------------------------- interface used for cfdemCloud -------------------------------------- */

  inline const std::vector<std::string>& forceModelList() const { return forceModelList_; }

  inline const std::vector<std::string>& momCoupleModelList() const { return momCoupleModelList_; }

  inline const std::vector<std::string>& liggghtsCommandModelList() const { return liggghtsCommandModelList_; }

  inline bool verbose() const { return verbose_; }

  inline bool solveFlow() const { return solveFlow_; }

  inline const std::string& modelType() const { return modelType_; }

  inline const std::string& turbulenceModelType() const { return turbulenceModelType_; }

  inline bool allowUseSubCFDTimeStep() const { return allowUseSubCFDTimeStep_; }

  inline int couplingInterval() const { return couplingInterval_; }

  inline bool checkPeriodicCells() const { return checkPeriodicCells_; }

  inline const Foam::vector& periodicCheckRange() const { return periodicCheckRange_; }

  inline const std::string& ddtVoidFractionType() const { return ddtVoidFractionType_; }

  /* --------------------------------- interface used for cfdemCloudIB ------------------------------------ */

  inline double refineMeshSkin() const { return refineMeshSkin_; }

  inline int refineMeshKeepInterval() const { return refineMeshKeepInterval_; }

  /* --------------------------------- interface used for cfdemCloudMix ----------------------------------- */

  inline double expandedCellScale() const { return expandedCellScale_; }

  inline bool checkFineParticle(double ratio) const { return ratio > 0 && ratio >= fineParticleRatio_; }

  inline bool checkMiddleParticle(double ratio) const {
    return ratio > 0 && (coarseParticleRatio_ <= ratio && ratio < fineParticleRatio_);
  }

  inline bool checkCoarseParticle(double ratio) const { return ratio > 0 && ratio < coarseParticleRatio_; }

 protected:
  const fvMesh& mesh_;

  //! \note 在当前类中一定要最先声明 couplingPropertiesDict_ 和 liggghtsCommandsDict_
  const IOdictionary& couplingPropertiesDict_;

  const IOdictionary& liggghtsCommandsDict_;

  //! \brief 颗粒尺度系数
  double fineParticleRatio_;

  double coarseParticleRatio_;

  //! \brief 颗粒扩展网格系数
  double expandedCellScale_;

  //! \brief dict 中指定的所有 force model 的名称
  std::vector<std::string> forceModelList_;

  //! \brief dict 中指定的所有 mom couple model 的名称
  std::vector<std::string> momCoupleModelList_;

  //! \brief dict 中指定的所有 liggghts command model 的名称
  std::vector<std::string> liggghtsCommandModelList_;

  //! \brief 是否打印多余调试信息
  bool verbose_;

  //! \brief 是否需要求解流体，在求解器代码中使用，如果为 false，则会跳过流体方程的计算
  bool solveFlow_;

  //! \brief 模型类型(A, B, Bfull, none)
  std::string modelType_;

  //! \brief dict 中指定的湍流模型的名称
  std::string turbulenceModelType_;

  //! \brief 是否允许耦合间隔中存在多个 CFD 时间步长（通常设置为 false，则 CFD 时间步长 == 耦合时间步长）
  bool allowUseSubCFDTimeStep_;

  //! \brief 耦合间隔，单位：DEM 时间步，couplingInterval_ * DEMts_ 就是耦合间隔的秒数
  int couplingInterval_;

  //! \brief used to de-activate mirroring across periodic boundary conditions.
  bool checkPeriodicCells_;

  //! \brief default = (1,1,1), periodic checks will be done according to periodicCheckRange_
  const Foam::vector periodicCheckRange_;

  /*!
   * \brief 重构网格的颗粒直径因数，在颗粒 refineMeshSkin_ * diameter 中的网格设置 interface 用于 refine meshs
   * \note used for cfdemCloudIB or cfdemCloudMix
   */
  double refineMeshSkin_;

  /*!
   * \brief 重构网格保留时间间隔(以流体时间步为单位)
   * \note 如果设置为 100, 即重构网格会保留 100 个时间步，默认值为 0
   */
  int refineMeshKeepInterval_;

  //! \brief ddtVoidFraction的类型，可以是"a" "b" or "off"，默认是"off"
  std::string ddtVoidFractionType_;

#if CFDEM_MIX_CLOUD

 protected:
  //! \brief 是否用于求解器 cfdemSolverIB
  bool usedForSolverIB_;

  //! \brief 是否用于求解器 cfdemSolverPiso
  bool usedForSolverPiso_;

  //! \brief 是否使用 dynamic refine mesh
  bool useDynamicRefineMesh_;

  //! \brief true if particle is fixed
  bool fixedParticle_;

  //! \brief 来流速度
  vector flowVelocity_;

 public:
  inline bool checkCoarseParticle(const double& ratio) const { return ratio > 0 && ratio < coarseParticleRatio_; }
  inline bool checkMiddleParticle(const double& ratio) const {
    return ratio > 0 && (coarseParticleRatio_ <= ratio && ratio < fineParticleRatio_);
  }
  inline bool checkFineParticle(const double& ratio) const { return ratio > 0 && ratio >= fineParticleRatio_; }
  inline bool checkFineAndMiddleParticle(const double& ratio) const { return !checkCoarseParticle(ratio); }
  /*!
   * \brief 是否需设置 field for coarse particle
   * \param index            <[in] 颗粒索引
   * \param dimensionRatios  <[in] 颗粒尺度比值
   */
  inline bool needSetFieldForCoarseParticle(const int& index, const std::vector<double> dimensionRatios) const {
    if (usedForSolverIB_ || useDynamicRefineMesh_) {
      return true;
    } else if (!dimensionRatios.empty() && !usedForSolverPiso_) {
      return index < static_cast<int>(dimensionRatios.size()) ? checkCoarseParticle(dimensionRatios[index]) : false;
    } else {
      return false;
    }
  }
  inline bool usedForSolverIB() const { return usedForSolverIB_; }
  inline bool usedForSolverPiso() const { return usedForSolverPiso_; }
  inline bool useDynamicRefineMesh() const { return useDynamicRefineMesh_; }
  inline bool fixedParticle() const { return fixedParticle_; }
  inline const vector& flowVelocity() const { return flowVelocity_; }
#endif  // CFDEM_MIX_CLOUD
};

}  // namespace Foam

#endif  // __COUPLING_PROPERTIES_H__
