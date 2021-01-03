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
\*---------------------------------------------------------------------------*/

#include "./mom_couple_model.h"

namespace Foam {

cfdemDefineTypeName(momCoupleModel);

cfdemDefineNewFunctionMap(momCoupleModel);

cfdemDefineConstructNewFunctionMap(momCoupleModel);

cfdemDefineDestroyNewFunctionMap(momCoupleModel);

cfdmeDefineBaseTypeNewWithTypeName(std::unique_ptr, momCoupleModel,
                                   (cfdemCloud & cloud, const dictionary& dict, const std::string& modelName),
                                   modelName, (cloud));

//! \brief Constructor
momCoupleModel::momCoupleModel(cfdemCloud& cloud) : cloud_(cloud), maxAlpha_(1 - Foam::SMALL) {}

//! \brief Destructor
momCoupleModel::~momCoupleModel() {}

//! \brief get implicit momentum source field
tmp<volScalarField> momCoupleModel::impMomSource() {
  FatalError << "the solver calls for impMomSource(), please set 'momCoupleModel' to type 'implicitCouple'"
             << abort(FatalError);
  tmp<volScalarField> source;
  return source;
}

//! \brief get explicit momentum source field
tmp<volVectorField> momCoupleModel::expMomSource() {
  FatalError << "the solver calls for expMomSource(), please set 'momCoupleModel' to type 'explicitCouple'"
             << abort(FatalError);
  tmp<volVectorField> source;
  return source;
}

}  // namespace Foam
