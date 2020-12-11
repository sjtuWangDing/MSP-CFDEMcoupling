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
  The IB voidFraction model is supposed to be used when a particle
  (or its representation) is bigger than a CFD cell. The voidfraction field
  is set in those cell whose centres are inside the particle. The model is
  specially designed for cfdemSolverIB and creates a smooth transition of
  the voidfraction at the particle surface. Cells which are only partially
  covered by solid are marked by voidfraction values between 0 and 1
  respectively.

  The region of influence of a particle can be increased artificially
  by “scaleUpVol”, which blows up the particles, but keeps their volume
  (for voidfraction calculation) constant.

Syntax
  voidfractionModel IB;
  IBProps
  {
    maxCellsPerParticle number1;
    alphaMin number2;
    scaleUpVol number3;
    checkPeriodicCells ;
  }
  - number1 = maximum number of cells covered by a particle (search will
              fail when more than number1 cells are covered by the particle)
  - number2 = minimum limit for voidfraction
  - number3 = diameter of the particle’s representation is artificially
              increased according to number3 * Vparticle, volume remains
              unaltered!
  - checkPeriodicCells = (optional, default false) flag for considering the
              minimal distance to all periodic images of this particle

Class
  IBVoidFraction
\*---------------------------------------------------------------------------*/

#ifndef __IB_VOID_FRACTION_H__
#define __IB_VOID_FRACTION_H__

#include "./void_fraction_model.h"

namespace Foam {

class IBVoidFraction : public voidFractionModel {

public:

  cfdemTypeName("IB")

  cfdemDefineNewFunctionAdder(voidFractionModel, IBVoidFraction)

  //! \brief Constructor
  IBVoidFraction(cfdemCloud& cloud);

  //! \brief Destructor
  virtual ~IBVoidFraction();

protected:

  /*!
   * \brief 设置单个颗粒的体积分数场
   * \param index 颗粒索引
   */
  void setVolumeFractionForSingleParticle(int index);

  /*!
   * \brief 构建颗粒覆盖的所有网格的哈希集合
   * \note 设置为递归函数,  通过哈希器将网格编号转换为哈希值,  并存入 set 中以便于搜索
   * \param cellID         <[in] 递归循环中要检索网格编号
   * \param particleCentre <[in] 颗粒中心位置
   * \param radius         <[in] 颗粒半径
   * \param hashSett       <[in, out] 需要构建的哈希集
   */
  void buildLabelHashSetForVolumeFraction(const label cellID,
                                          const Foam::vector& particleCentre,
                                          const double radius,
                                          labelHashSet& hashSett);

  /*!
   * \brief 计算距离系数，对任意一个网格, 如果网格中心 c 在颗粒内部, 但是它的某个角点 p
   *   不在颗粒内部, 则计算 c 与 p 的连线与颗粒表面的交点 i 到网格中心 c 的距离, 即
   *   求解 x 的二元一次方程
   *   (x * (vector_p - vector_c) - vector_particle) & 
   *   (x * (vector_p - vector_c) - vector_particle) == radius * radius
   *   等价于函数体中定义的: a*(x^2) - b*x + c = 0
   * \param radius         <[in] 颗粒半径
   * \param particleCentre <[in] 颗粒中心
   * \param pointInside    <[in] 网格中心
   * \param pointOutside   <[in] 网格角点
   */
  double segmentParticleIntersection(double radius,
                                     const Foam::vector& particleCentre,
                                     const Foam::vector& pointInside,
                                     const Foam::vector& pointOutside) const;

  /*!
   * \brief 获取 Corona Point
   * \note Corona Point 是在以网格中心为中心, 半径为 corona 的球面上, 距离颗粒中心最远的点
   *       其中, 半径 corona = 0.5 * sqrt(3) * 网格等效半径
   *       网格等效半径 = pow(cellVolume, 1.0 / 3.0)
   * \param particleCentre  <[in] 指定颗粒中心
   * \param cellCentre      <[in] 指定网格中心
   * \param corona          <[in] 指定网格的等效半径
   */
  Foam::vector getCoronaPointPosition(const Foam::vector& particleCentre,
                                      const Foam::vector& cellCentre,
                                      const scalar corona) const;

private:
  dictionary subPropsDict_;

  const double alphaMin_;

  const double alphaMax_;
};

} // namespace Foam

#endif // __IB_VOID_FRACTION_H__
