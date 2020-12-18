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
  cfdemCloudIB derived from cfdemCloud

Class
  Foam::cfdemCloudIB
\*---------------------------------------------------------------------------*/

#ifndef __CFDEM_CLOUD_IB_H__
#define __CFDEM_CLOUD_IB_H__

#include "cloud/cfdem_cloud.h"

namespace Foam {

class cfdemCloudIB : public cfdemCloud {
 public:
  //! \brief Constructed from mesh
  cfdemCloudIB(const fvMesh& mesh);

  //! \brief Destructor
  ~cfdemCloudIB();

  /*!
   * \brief 更新函数
   * \note used for cfdemSolverIB
   * \param volumeFraction  <[in, out] 大颗粒体积分数
   * \param interFace       <[in, out] 界面场，用于 dynamic mesh
   */
  void evolve(volScalarField& volumeFraction, volScalarField& interface);

  //! \brief 重新分配内存
  void reallocate();

  //! \brief 重新分配内存
  void reallocateTensorVector();

  //! \brief 从 DEM 获取数据
  void getDEMData();

  //! \brief 传递数据到 DEM
  void giveDEMData() const;

  void getDimensionRatio();

  /*!
   * \brief 更新网格，如果 mesh 是 Foam::dynamicRefineFvMesh 类型，则更新网格，
   *   如果是 Foam::staticFvMesh 或者其他类型，则不更新
   */
  void updateMesh(volScalarField& interface);

  //! @brief 确定颗粒周围 refined 网格的区域
  void setInterface(volScalarField& interface, const double scale = 1.0) const;

  //! @brief 确定颗粒周围 refined 网格的区域(每个方向的尺寸都是颗粒尺寸的两倍)
  void setInterface(volScalarField& interface, volScalarField& refineMeshKeepStep) const;

 public:
  inline double refineMeshSkin() const { return cProps_.refineMeshSkin(); }

  inline int refineMeshKeepInterval() const { return cProps_.refineMeshKeepInterval(); }

  inline bool meshHasUpdated() const { return meshHasUpdated_; }

  inline void setMeshHasUpdated(bool meshHasUpdated) { meshHasUpdated_ = meshHasUpdated; }

 protected:
  static const double particleMeshScale_;
  /*!
   * \brief 判断 mesh 是否被更新过
   * \note 在求解器中使用 dynamic mesh，如果 mesh 更新，则 mesh.update() 返回 true，否则返回 false
   */
  bool meshHasUpdated_;
};

const double cfdemCloudIB::particleMeshScale_ = 1.5;

}  // namespace Foam

#endif  // __CFDEM_CLOUD_IB_H__
