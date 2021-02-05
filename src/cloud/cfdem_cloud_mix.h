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
  cfdemCloudMix derived from cfdemCloud

Class
  Foam::cfdemCloudMix
\*---------------------------------------------------------------------------*/

#ifndef __CFDEM_CLOUD_MIX_H__
#define __CFDEM_CLOUD_MIX_H__

#include "cloud/cfdem_cloud.h"

namespace Foam {

class cfdemCloudMix : public cfdemCloud {
 public:
  //! \brief Constructed from mesh
  cfdemCloudMix(const fvMesh& mesh);

  //! \brief Destructor
  ~cfdemCloudMix();

  //! \brief Runtime type information
  cfdemTypeName("cfdemCloudMix");

  //! \brief 重新分配内存
  void reallocate();

  //! \brief 从 DEM 获取数据
  void getDEMData();

  void printParticleInfo() const;

  double expandedCellScale() const { return cProps_.expandedCellScale(); }

  bool checkFineParticle(int index) const { return cProps_.checkFineParticle(getDimensionRatio(index)); }

  bool checkMiddleParticle(int index) const { return cProps_.checkMiddleParticle(getDimensionRatio(index)); }

  bool checkCoarseParticle(int index) const { return cProps_.checkCoarseParticle(getDimensionRatio(index)); }

  /*!
   * \brief 更新函数
   * \param U          <[in] 流体速度场
   * \param voidF      <[in, out] 小颗粒空隙率场
   * \param volumeF    <[in, out] 大颗粒空隙率场
   * \param Us         <[in, out] 局部平均速度场
   * \param Ksl        <[in, out] 动量交换场
   * \param interface  <[in, out] 界面场
   */
  void evolve(volVectorField& U, volScalarField& voidF, volScalarField& volumeF, volVectorField& Us,
              volScalarField& Ksl, volScalarField& interface);
};

}  // namespace Foam

#endif  // __CFDEM_CLOUD_MIX_H__
