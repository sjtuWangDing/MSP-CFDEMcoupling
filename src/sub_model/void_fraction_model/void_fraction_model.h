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
  The voidFractionModel is the base class for models to represent the DEM
  particle’s volume in the CFD domain via a voidFraction field.

Syntax
  Defined in couplingProperties dictionary:
  voidfractionModel model;

Class
  voidFractionModel, default is noVoidFraction
\*---------------------------------------------------------------------------*/

#ifndef __VOID_FRACTION_MODEL_H__
#define __VOID_FRACTION_MODEL_H__

#include "base/run_time_selection_tables.h"
#include "cloud/cfdem_cloud.h"

namespace Foam {

//! \brief 空隙率模型
class voidFractionModel {
 public:
  //! \brief Runtime type information
  cfdemTypeName("voidFractionModel");

  //! \brief Declare runtime constructor selection
  cfdemDeclareRunTimeSelection(autoPtr, voidFractionModel, (cfdemCloud & cloud), (cloud));

  //! \brief Selector
  static autoPtr<voidFractionModel> New(cfdemCloud& cloud, const dictionary& dict);

  //! \brief Constructor
  voidFractionModel(cfdemCloud& cloud);

  //! \brief Destructor
  virtual ~voidFractionModel();

  //! \brief 计算颗粒尺寸与其周围网格平均尺寸的比值, 并将颗粒索引按照颗粒尺寸归类
  void getDimensionRatios(const base::CITensor1& findCellIDs, const base::CDTensor1& dimensionRatios) const;

  //! \brief 重置空隙率
  inline void resetVoidFraction() {
    voidFractionPrev_ == voidFractionNext_;
    voidFractionNext_ == dimensionedScalar("one", voidFractionNext_.dimensions(), 1.);
  }

  //! \brief 重置体积分数
  inline void resetVolumeFraction() {
    volumeFractionPrev_ == volumeFractionNext_;
    volumeFractionNext_ == dimensionedScalar("one", volumeFractionNext_.dimensions(), 1.);
  }

  /*!
   * \brief 判断某一点是否在颗粒中
   * \param point
   * \param particlePos 颗粒中心坐标
   * \param radius 颗粒半径
   * \param scale 半径因子
   * \return < 0: 在颗粒中，范围为 [-1, 0)
   *         > 0: 在颗粒外
   *         == 0: 在颗粒表面
   */
  inline double pointInParticle(const Foam::vector& point, const Foam::vector& particlePos, const double radius,
                                const double scale = 1.0) const {
    return Foam::magSqr(point - particlePos) / (scale * scale * radius * radius) - 1.0;
  }

  inline const volScalarField& voidFractionPrev() const { return voidFractionPrev_; }

  inline const volScalarField& voidFractionNext() const { return voidFractionNext_; }

  inline const volScalarField& volumeFractionPrev() const { return volumeFractionPrev_; }

  inline const volScalarField& volumeFractionNext() const { return volumeFractionNext_; }

 protected:
  /*!
   * \brief 构建颗粒覆盖的所有网格的哈希集合
   * \note 设置为递归函数，通过哈希器将网格编号转换为哈希值，并存入 set 中以便于搜索
   * \param hashSett    <[in, out] 需要构建的哈希集
   * \param cellID      <[in] 递归循环中要检索网格编号
   * \param particlePos <[in] 颗粒中心位置
   * \param radius      <[in] 颗粒半径
   * \param scale       <[in] 颗粒半径扩大系数
   */
  void buildLabelHashSetForCoveredCell(labelHashSet& hashSett, const int cellID, const Foam::vector& particlePos,
                                       const double radius, const double scale) const;

 protected:
  cfdemCloud& cloud_;

  //! \brief 小颗粒空隙率
  volScalarField voidFractionPrev_;

  volScalarField voidFractionNext_;

  //! \brief 大颗粒体积分数
  volScalarField volumeFractionPrev_;

  volScalarField volumeFractionNext_;

  //! \brief fine颗粒覆盖的网格数量的上限
  int maxCellsNumPerFineParticle_;

  //! \brief middle 颗粒覆盖的网格数量的上限
  int maxCellsNumPerMiddleParticle_;

  //! \brief coarse 颗粒覆盖的网格数量的上限
  int maxCellsNumPerCoarseParticle_;
};

}  // namespace Foam

#endif  // __VOID_FRACTION_MODEL_H__
