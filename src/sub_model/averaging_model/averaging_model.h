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
  The averaging model performs the Lagrangian->Eulerian mapping of
  data (e.g. particle velocities).

Syntax
  averagingModel modelName;

Class
  averagingModel
\*---------------------------------------------------------------------------*/

#ifndef __AVERAGING_MODEL_H__
#define __AVERAGING_MODEL_H__

#include "base/run_time_selection_tables.h"
#include "cloud/cfdem_cloud.h"

namespace Foam {

//! \brief averagingModel
class averagingModel {
 public:
  //! \brief Runtime type information
  cfdemTypeName("averagingModel");

  //! \brief Declare runtime constructor selection
  cfdemDeclareRunTimeSelection(autoPtr, averagingModel, (cfdemCloud & cloud), (cloud));

  //! \brief Selector
  static autoPtr<averagingModel> New(cfdemCloud& cloud, const dictionary& dict);

  //! \brief Constructor
  averagingModel(cfdemCloud& cloud);

  //! \brief Destructor
  virtual ~averagingModel();

  inline void resetUs() {
    UsPrev_ == UsNext_;
    UsNext_ == dimensionedVector("zero", UsNext_.dimensions(), vector::zero);
  }

  inline void resetUsWeightField() { UsWeightField_ == dimensionedScalar("zero", UsWeightField_.dimensions(), 0.0); }

  inline const volVectorField& UsPrev() const { return UsPrev_; }

  inline const volVectorField& UsNext() const { return UsNext_; }

  inline const volScalarField& UsWeightField() const { return UsWeightField_; }

 protected:
  cfdemCloud& cloud_;

  //! @brief 上个时间步中局部平均颗粒速度场
  volVectorField UsPrev_;

  //! @brief 下个时间步中局部平均颗粒速度场
  volVectorField UsNext_;

  //! @brief 覆盖网格的所有颗粒对网格的影响系数总和
  volScalarField UsWeightField_;
};

}  // namespace Foam

#endif  // __AVERAGING_MODEL_H__
