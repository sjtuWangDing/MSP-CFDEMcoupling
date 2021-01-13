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
  The momCoupleModel is the base class for momentum exchange between
  DEM and CFD simulation.

Syntax
  Defined in couplingProperties dictionary:
  momCoupleModel
  (
    model
  );

Class
  Foam::momCoupleModel
\*---------------------------------------------------------------------------*/

#ifndef __MOM_COUPLE_MODEL_H__
#define __MOM_COUPLE_MODEL_H__

#include "base/run_time_selection_tables.h"
#include "cloud/cfdem_cloud.h"

namespace Foam {

class momCoupleModel {
 public:
  //! \brief Runtime type information
  cfdemBaseTypeName("momCoupleModel", "");

  //! \brief Declare runtime constructor selection
  cfdemDeclareRunTimeSelection(std::unique_ptr, momCoupleModel, (cfdemCloud & cloud), (cloud));

  //! \brief Selector
  static std::unique_ptr<momCoupleModel> New(cfdemCloud& cloud, const dictionary& dict, const std::string& modelName);

  //! \brief Constructor
  momCoupleModel(cfdemCloud& cloud);

  //! \brief Destructor
  virtual ~momCoupleModel();

  //! \brief reset mom field
  virtual void resetMomSourceField() = 0;

  //! \brief get implicit momentum source field
  virtual tmp<volScalarField> impMomSource();

  //! \brief get explicit momentum source field
  virtual tmp<volVectorField> expMomSource();

 protected:
  cfdemCloud& cloud_;

  //! \brief max fluid volume fraction to calculate exchange field
  double maxAlpha_;
};

}  // namespace Foam

#endif  // __MOM_COUPLE_MODEL_H__
