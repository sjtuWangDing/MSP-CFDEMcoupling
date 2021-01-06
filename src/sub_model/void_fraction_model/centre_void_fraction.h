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
  The centre voidFraction model calculates the voidfraction in a CFD cell
  accounting for the volume of the particles whose centres are inside the
  cell. The particle volume occupied in the CFD domain can be adjusted
  by the parameter “weight”, using Vparticle = dsphere^3*pi/6*weight.

Syntax
  voidfractionModel centre;
  centreProps
  {
    alphaMin number1;
    weight number2;
  }
  - number1 = minimum limit for voidfraction
  - number2 = (optional) scaling of the particle volume to account for
              porosity or agglomerations.

Class
  noVoidFraction
\*---------------------------------------------------------------------------*/

#ifndef __CENTRE_VOID_FRACTION_H__
#define __CENTRE_VOID_FRACTION_H__

#include "./void_fraction_model.h"

namespace Foam {

class centreVoidFraction : public voidFractionModel {
 public:
  cfdemTypeName("centre");

  cfdemDefineNewFunctionAdder(voidFractionModel, centreVoidFraction);

  //! \brief Constructor
  centreVoidFraction(cfdemCloud& cloud);

  //! \brief Destructor
  ~centreVoidFraction();

  void setVoidFraction();

  //! \brief 输出空隙率相关信息
  void printVoidFractionInfo() const {}

 protected:
  //! \brief 设置索引为 index 的单个颗粒的空隙率
  void setVoidFractionForSingleParticle(const int index, std::unordered_map<int, Foam::vector>& parMap);

 private:
  dictionary subPropsDict_;

  //! \brief 空隙率的最小值
  double alphaMin_;
};

}  // namespace Foam

#endif  // __CENTRE_VOID_FRACTION_H__
