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
  Foam::noAveragingModel

Syntax
  averagingModel noAveragingModel;
\*---------------------------------------------------------------------------*/

#ifndef __NO_AVERAGING_MODEL_H__
#define __NO_AVERAGING_MODEL_H__

#include "./averaging_model.h"

namespace Foam {

//! \brief noAveragingModel
class noAveragingModel : public averagingModel {
 public:
  //! \brief Runtime type information
  cfdemTypeName("noAveragingModel");

  cfdemDefineNewFunctionAdder(averagingModel, noAveragingModel);

  //! \brief Constructor
  noAveragingModel(cfdemCloud& cloud);

  //! \brief Destructor
  ~noAveragingModel();

  /*!
   * \brief 计算局部平均矢量场
   * \param valueField   <[in, out] 需要被局部平均化场
   * \param weightField  <[in, out] 权重系数平均化场
   * \param value        <[in] 用于局部平均化的颗粒数据(lagrange value)
   * \param weight       <[in] 用于局部平均化的权重系数(lagrange value)
   */
  void setVectorFieldAverage(volVectorField& valueField, volScalarField& weightField, const base::CDExTensor2& value,
                             const std::vector<base::CDTensor1>& weight) {}
};

}  // namespace Foam

#endif  // __NO_AVERAGING_MODEL_H__
