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
  The divided voidFraction model is supposed to be used when a particle
  (or its representation) is in the size range of a CFD cell. The particle
  has radius R and it’s volume is divided in 29 non-overlapping regions of
  equal volume.

Syntax
  voidFractionModel divided;
  dividedProps
  {
    alphaMin number1;
    interpolation;
    weight number2;
    porosity number3;
    procBoundaryCorrection Switch1;
    verbose;
  }
  - number1 = minimum limit for voidFraction
  - interpolation = flag to interpolate voidFraction to particle
                    positions (normally off)
  - number2 = (optional) scaling of the particle volume to account
              for porosity or agglomerations.
  - number3 = (optional) diameter of the particle’s representation
              is artificially increased according to number2 * Vparticle,
              volume remains unaltered!
  - Switch1 = (optional, default false) allow for correction at
              processor boundaries. This requires the use of engineIB
              and vice versa.
  - verbose = (optional, default false) flag for debugging output

Class
  Foam::divided
\*---------------------------------------------------------------------------*/

#ifndef __DIVIDED_VOID_FRACTION_H__
#define __DIVIDED_VOID_FRACTION_H__

#include "./void_fraction_model.h"

namespace Foam {

class dividedVoidFraction : public voidFractionModel {
 public:
  cfdemTypeName("divided");

  cfdemDefineNewFunctionAdder(voidFractionModel, dividedVoidFraction);

  //! \brief Constructor
  dividedVoidFraction(cfdemCloud& cloud);

  //! \brief Destructor
  virtual ~dividedVoidFraction();

  /*!
   * \brief 当颗粒的体积在CFD单元的尺寸范围内时，应使用分割(divided)空隙率模型
   * \note 粒子的半径为R，颗粒的体积被划分为29个等体积的不重叠区域
   */
  static const int numberOfMarkerPoints_ = 29;

  //! \brief 计算空隙率
  void setVoidFraction();

  //! \brief 输出空隙率相关信息
  void printVoidFractionInfo() const;

 protected:
  //! \brief 设置索引为 index 的单个颗粒的空隙率
  void setVoidFractionForSingleParticle(const int index, std::unordered_map<int, Foam::vector>& parMap);

 private:
  dictionary subPropsDict_;

  //! \brief 是否打印详细信息
  bool verbose_;

  //! \brief 空隙率的最小值
  double alphaMin_;

  //! \brief 由于空隙率限制而导致颗粒体积的累积损失
  double tooMuch_;

  /*!
   * \brief 标志点坐标
   * \note 这里 offsets_ = (标志点中心到颗粒中心的矢量 / 半径)，不是实际坐标
   */
  Foam::vector offsets_[numberOfMarkerPoints_];
};

}  // namespace Foam

#endif  // __DIVIDED_VOID_FRACTION_H__
