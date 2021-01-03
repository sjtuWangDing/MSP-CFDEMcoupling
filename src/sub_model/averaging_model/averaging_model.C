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
  averagingModel
\*---------------------------------------------------------------------------*/

#include "./averaging_model.h"

namespace Foam {

cfdemDefineTypeName(averagingModel);

cfdemDefineNewFunctionMap(averagingModel);

cfdemDefineConstructNewFunctionMap(averagingModel);

cfdemDefineDestroyNewFunctionMap(averagingModel);

cfdmeDefineBaseTypeNew(autoPtr, averagingModel, (cfdemCloud & cloud, const dictionary& dict), dict, (cloud));

//! \brief Constructor
averagingModel::averagingModel(cfdemCloud& cloud)
    : cloud_(cloud),
      UsPrev_(IOobject("UsPrev", cloud.mesh().time().timeName(), cloud.mesh(), IOobject::READ_IF_PRESENT,
                       IOobject::AUTO_WRITE),
              cloud_.mesh(), dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0), Foam::vector::zero)),
      UsNext_(IOobject("UsNext", cloud.mesh().time().timeName(), cloud.mesh(), IOobject::READ_IF_PRESENT,
                       IOobject::AUTO_WRITE),
              cloud_.mesh(), dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0), Foam::vector::zero)),
      UsWeightField_(IOobject("UsWeightField", cloud.mesh().time().timeName(), cloud.mesh(), IOobject::NO_READ,
                              IOobject::AUTO_WRITE),
                     cloud.mesh(), dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)) {}

//! \brief Destructor
averagingModel::~averagingModel() {}

}  // namespace Foam
