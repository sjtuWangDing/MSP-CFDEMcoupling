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
  ShirgaonkarIB
\*---------------------------------------------------------------------------*/

#include "./Shirgaonkar_IB.h"

namespace Foam {

cfdemDefineTypeName(ShirgaonkarIB)

cfdemCreateNewFunctionAdder(forceModel, ShirgaonkarIB)

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
ShirgaonkarIB::ShirgaonkarIB(cfdemCloud& cloud)
  : forceModel(cloud),
    subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
    velFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("velFieldName", "U").c_str()),
    pressureFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("pressureFieldName", "p").c_str()),
    U_(cloud.mesh().lookupObject<volVectorField>(velFieldName_)),
    p_(cloud.mesh().lookupObject<volScalarField>(pressureFieldName_)),
    useTorque_(subPropsDict_.lookupOrDefault<bool>("useTorque", false)) {
  createForceSubModels(subPropsDict_, kResolved);
}

ShirgaonkarIB::~ShirgaonkarIB() {}

void ShirgaonkarIB::setForce() {
  Info << "Setting ShirgaonkarIB force..." << endl;
  volVectorField IBDrag = forceSubModel_->IBDrag(U_, p_);
  Foam::vector drag = Foam::vector::zero;
  Foam::vector torque = Foam::vector::zero;
  Foam::vector cellPos = Foam::vector::zero;

  for (int index = 0; index <= cloud_.numberOfParticles(); ++index) {
    // init
    drag = Foam::vector::zero;
    torque = Foam::vector::zero;
    // get index's particle center position
    Foam::vector particleCenterPos = cloud_.getPosition(index);
    // loop all mesh of current particle
    for (int subCell = 0; subCell < cloud_.particleOverMeshNumber()[index]; ++subCell) {
      label cellI = cloud_.cellIDs()[index][subCell];
      if (cellI > -1) { // cell Found
        cellPos = cloud_.mesh().C()[cellI];
        drag += IBDrag[cellI] * IBDrag.mesh().V()[cellI];
        torque += (cellPos - particleCenterPos) ^ IBDrag[cellI] * IBDrag.mesh().V()[cellI];
      }
    }
    // write particle data to global array
    // index - particle index
    // drag - total drag
    forceSubModel_->partToArray(index, drag, Foam::vector::zero, Foam::vector::zero, 0);

    if (forceSubModel_->verbose()) {
      Info << "drag on particle " << index << ": ["
        << drag[0] << ", " << drag[1] << ", " << drag[2] << "]" << endl;
    }

    if (useTorque_) {
      forceSubModel_->addTorque(index, torque);
    }
  }
}

} // namespace Foam
