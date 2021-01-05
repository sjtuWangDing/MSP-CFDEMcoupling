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
  This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
  and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

Class
  forceModel
\*---------------------------------------------------------------------------*/

#ifndef __FORCE_MODEL_H__
#define __FORCE_MODEL_H__

#include "./force_sub_model.h"
#include "base/run_time_selection_tables.h"
#include "cloud/cfdem_cloud.h"
#include "interpolationCellPointFace.H"

namespace Foam {

//! \brief forceModel
class forceModel {
 public:
  //! \brief Runtime type information
  cfdemTypeName("forceModel");

  //! \brief Declare runtime constructor selection
  cfdemDeclareRunTimeSelection(std::unique_ptr, forceModel, (cfdemCloud & cloud), (cloud));

  //! \brief Selector
  static std::unique_ptr<forceModel> New(cfdemCloud& cloud, const dictionary& dict, const std::string& modelName);

  //! \brief Constructor
  forceModel(cfdemCloud& cloud);

  //! \brief Destructor
  virtual ~forceModel();

  virtual void setForce() {
    FatalError << "forceModel:setForce(): using base class function, please use derived class function\n"
               << abort(FatalError);
  }

  /*!
   * \brief create forceSubModel_
   * \param subPropsDict the dictionary of current force model
   * \param forceType force type
   */
  void createForceSubModels(const dictionary& subPropsDict, EForceType forceType);

  inline std::shared_ptr<forceSubModel> forceSubM() const { return forceSubModel_; }

 protected:
  cfdemCloud& cloud_;

  //! \brief 当前 force model 的 forceSubModel
  std::shared_ptr<forceSubModel> forceSubModel_;

  //! \brief 是否激活探针
  bool useProbe_;

  autoPtr<interpolation<Foam::vector>> UInterpolator_;

  autoPtr<interpolation<Foam::scalar>> voidFractionInterpolator_;
};

}  // namespace Foam

#endif  // __FORCE_MODEL_H__
