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

Class
  Foam::mix
\*---------------------------------------------------------------------------*/

#ifndef __MIX_VOID_FRACTION_H__
#define __MIX_VOID_FRACTION_H__

#include "./void_fraction_model.h"

namespace Foam {

class mixVoidFraction : public voidFractionModel {
 public:
  cfdemTypeName("mix");

  cfdemDefineNewFunctionAdder(voidFractionModel, mixVoidFraction);

  //! \brief Constructor
  mixVoidFraction(cfdemCloud& cloud);

  //! \brief Destructor
  virtual ~mixVoidFraction();

  //! \brief 计算空隙率
  void setVoidFraction();

  //! \brief 输出空隙率相关信息
  void printVoidFractionInfo() const;

 protected:
  //! \brief 设置索引为 index 的单个颗粒的空隙率
  //! \note used for fine particle
  void setVoidFractionForSingleParticle(const int index, const int findCellID,
                                        std::unordered_map<int, Foam::vector>& parMap);

  /*!
   * \brief 设置单个颗粒的体积分数场
   * \param index 颗粒索引
   * \param set   颗粒覆盖的网格索引的集合
   */
  void setVolumeFractionForSingleParticle(int index, std::unordered_set<int>& set);

  /*!
   * \brief 构建颗粒覆盖的所有网格的哈希集合
   * \note 设置为递归函数,  通过哈希器将网格编号转换为哈希值,  并存入 set 中以便于搜索
   * \param cellID         <[in] 递归循环中要检索网格编号
   * \param particleCentre <[in] 颗粒中心位置
   * \param radius         <[in] 颗粒半径
   * \param set            <[in, out] 颗粒覆盖的网格索引的集合
   */
  void buildSetForVolumeFraction(const label cellID, const Foam::vector& particleCentre, const double radius,
                                 std::unordered_set<int>& set);

 private:
  dictionary subPropsDict_;

  //! \brief 是否打印详细信息
  bool verbose_;

  //! \brief 空隙率的最小值
  double alphaMin_;

  //! \brief 由于空隙率限制而导致颗粒体积的累积损失
  double tooMuch_;
};

}  // namespace Foam

#endif  // __MIX_VOID_FRACTION_H__
