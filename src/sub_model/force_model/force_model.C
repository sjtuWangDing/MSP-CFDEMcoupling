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
  demoModel
\*---------------------------------------------------------------------------*/

#include "./force_model.h"

namespace Foam {

cfdemDefineTypeName(forceModel)

    cfdemDefineNewFunctionMap(forceModel)

        cfdemDefineConstructNewFunctionMap(forceModel)

            cfdemDefineDestroyNewFunctionMap(forceModel)

                cfdmeDefineBaseTypeNewWithTypeName(std::unique_ptr, forceModel,
                                                   (cfdemCloud & cloud, const dictionary& dict,
                                                    const std::string& modelName),
                                                   modelName, (cloud))

    //! \brief Constructor
    forceModel::forceModel(cfdemCloud& cloud)
    : cloud_(cloud),
      forceSubModel_(nullptr),
      useProbe_(false),
      impParticleForces_(IOobject("impParticleForces", cloud.mesh().time().timeName(), cloud.mesh(),
                                  IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                         cloud.mesh(), dimensionedVector("zero", dimensionSet(1, 1, -2, 0, 0),
                                                         vector(0, 0, 0))  // [N] == [kg * m / s^2]
                         ),
      expParticleForces_(IOobject("expParticleForces", cloud.mesh().time().timeName(), cloud.mesh(),
                                  IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                         cloud.mesh(), dimensionedVector("zero", dimensionSet(1, 1, -2, 0, 0),
                                                         vector(0, 0, 0))  // [N] == [kg * m / s^2]
                         ) {}

//! \brief Destructor
forceModel::~forceModel() {}

/*!
 * \brief create forceSubModel_
 * \param subPropsDict the dictionary of current force model
 * \param forceType force type
 */
void forceModel::createForceSubModels(const dictionary& subPropsDict, EForceType forceType) {
  // create force sub model
  forceSubModel_ = std::make_shared<forceSubModel>(cloud_, *this, subPropsDict);
  // read switches in force sub model
  forceSubModel_->readSwitches();
  // check switches
  forceSubModel_->checkSwitches(forceType);
}

}  // namespace Foam
