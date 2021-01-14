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

  // //! \brief 设置索引为 index 的单个颗粒的空隙率
  // //! \note used for middle particle
  // void setVoidFractionForSingleMiddleParticle(const int index, const int findMpiCellID,
  //                                             std::unordered_map<int, Foam::vector>& parMap);

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
