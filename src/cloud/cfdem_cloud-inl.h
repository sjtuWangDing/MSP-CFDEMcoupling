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
  cfdemCloud class managing DEM data for CFD-DEM coupling

Class
  Foam::cfdemCloud
\*---------------------------------------------------------------------------*/

#ifndef __CFDEM_CLOUD_INL_H__
#define __CFDEM_CLOUD_INL_H__

#include <memory>
#include <vector>

#if defined(version21) || defined(version16ext)
#include "turbulenceModel.H"
#elif defined(version15)
#include "RASModel.H"
#endif

// Need include header files of class which uses Foam::autoPtr to
// wrap model For liggghtsCommandModel, forceModel and momCoupleModel,
// not need to include header because they use share_ptr
#include "sub_model/averaging_model/averaging_model.h"
#include "sub_model/data_exchange_model/data_exchange_model.h"
#include "sub_model/locate_model/locate_model.h"
#include "sub_model/void_fraction_model/void_fraction_model.h"

namespace Foam {

inline const std::vector<std::shared_ptr<liggghtsCommandModel>>& cfdemCloud::liggghtsCommandModels() const {
  return liggghtsCommandModels_;
}

inline const std::vector<std::shared_ptr<forceModel>>& cfdemCloud::forceModels() const {
  return forceModels_;
}

inline const std::vector<std::shared_ptr<momCoupleModel>>& cfdemCloud::momCoupleModels() const {
  return momCoupleModels_;
}

inline std::shared_ptr<forceModel> cfdemCloud::forceM(int index) const {
  if (index < 0 || index >= static_cast<int>(forceModels_.size())) {
    FatalError << "forceM(): "
               << "index " << index << " out of boundary" << abort(FatalError);
  }
  return forceModels_[index];
}

// Foam::autoPtr<T> 中定义了 inline operator const T&() const;
inline const dataExchangeModel& cfdemCloud::dataExchangeM() const {
  return dataExchangeModel_;
}

inline const voidFractionModel& cfdemCloud::voidFractionM() const {
  return voidFractionModel_;
}

inline const locateModel& cfdemCloud::locateM() const {
  return locateModel_;
}

inline dataExchangeModel& cfdemCloud::dataExchangeM() {
  return dataExchangeModel_();
}

inline voidFractionModel& cfdemCloud::voidFractionM() {
  return voidFractionModel_();
}

inline locateModel& cfdemCloud::locateM() {
  return locateModel_();
}

#if defined(version24Dev)
inline const turbulenceModel& cfdemCloud::turbulence() const
#elif defined(version21) || defined(version16ext)
#ifdef compre
inline const compressible::turbulenceModel& cfdemCloud::turbulence() const
#else
inline const incompressible::turbulenceModel& cfdemCloud::turbulence() const
#endif
#elif defined(version15)
inline const incompressible::RASModel& cfdemCloud::turbulence() const
#endif
{
  return turbulence_;
}

}  // namespace Foam

#endif  // __CFDEM_CLOUD_INL_H__
