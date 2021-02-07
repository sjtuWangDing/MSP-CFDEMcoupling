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
  cfdemCloudIBOpti derived from cfdemCloud

Class
  Foam::cfdemCloudIBOpti
\*---------------------------------------------------------------------------*/

#ifndef __CFDEM_CLOUD_IB_OPTI_H__
#define __CFDEM_CLOUD_IB_OPTI_H__

#include "cloud/cfdem_cloud.h"

namespace Foam {

class cfdemCloudIBOpti : public cfdemCloud {
 public:
  //! \brief Constructed from mesh
  cfdemCloudIBOpti(const fvMesh& mesh);

  //! \brief Destructor
  ~cfdemCloudIBOpti();

  /*!
   * \brief 更新函数
   * \note used for cfdemSolverIB
   * \param volumeFraction  <[in, out] 大颗粒体积分数
   * \param interface       <[in, out] 界面场，用于 dynamic mesh
   */
  void evolve(volScalarField& volumeFraction, volScalarField& interface);

  void calcVelocityCorrection(volScalarField& p, volVectorField& U, volScalarField& phiIB,
                              volScalarField& voidfraction);

 protected:
  //! \brief 重新分配内存
  void reallocate();

  //! \brief 从 DEM 获取数据
  void getDEMData();

  //! \brief 传递数据到 DEM
  void giveDEMData() const;

  void printParticleInfo() const;

  /*!
   * \brief 更新网格，如果 mesh 是 Foam::dynamicRefineFvMesh 类型，则更新网格，
   *   如果是 Foam::staticFvMesh 或者其他类型，则不更新
   */
  void updateMesh(volScalarField& interface);

  //! \brief 确定颗粒周围 refined 网格的区域
  void setInterface(volScalarField& interface, const double scale = cfdemCloudIBOpti::particleMeshScale_) const;

 protected:
  static const double particleMeshScale_;

  int pRefCell_;

  double pRefValue_;
};

const double cfdemCloudIBOpti::particleMeshScale_ = 2.0;

}  // namespace Foam

#endif  // __CFDEM_CLOUD_IB_OPTI_H__
