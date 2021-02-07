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

cfdemDefineTypeName(ShirgaonkarIB);

cfdemCreateNewFunctionAdder(forceModel, ShirgaonkarIB);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
ShirgaonkarIB::ShirgaonkarIB(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      U_(cloud.globalF().U()),
      p_(cloud.globalF().p()),
      useTorque_(subPropsDict_.lookupOrDefault<bool>("useTorque", false)) {
  createForceSubModels(subPropsDict_, kResolved);
}

ShirgaonkarIB::~ShirgaonkarIB() {}

void ShirgaonkarIB::setForce() {
  Info << "Setting ShirgaonkarIB force..." << endl;
  volVectorField IBDrag = forceSubModel_->IBDrag(U_, p_);
  Foam::vector particleCenterPos = Foam::vector::zero;
  Foam::vector cellPos = Foam::vector::zero;
  Foam::vector drag = Foam::vector::zero;
  Foam::vector torque = Foam::vector::zero;
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkCoarseParticle(index)) {
      // init
      drag = Foam::vector::zero;
      torque = Foam::vector::zero;
      // get index's particle center position
      particleCenterPos = cloud_.getPosition(index);
      // loop all mesh of current particle
      for (int subCell = 0; subCell < cloud_.particleOverMeshNumber()[index]; ++subCell) {
        int cellID = cloud_.cellIDs()[index][subCell];
        if (cellID >= 0) {  // cell Found
          cellPos = cloud_.mesh().C()[cellID];
          drag += IBDrag[cellID] * IBDrag.mesh().V()[cellID];
          torque += (cellPos - particleCenterPos) ^ IBDrag[cellID] * IBDrag.mesh().V()[cellID];
        }
      }
      // write particle data to global array
      forceSubModel_->partToArray(index, drag, Foam::vector::zero, Foam::vector::zero, 0);

      if (forceSubModel_->verbose()) {
        Pout << "drag on particle " << index << ": [" << drag[0] << ", " << drag[1] << ", " << drag[2] << "]" << endl;
      }

      if (useTorque_) {
        forceSubModel_->addTorque(index, torque);
      }
    }
  }
  Info << "Setting ShirgaonkarIB - done" << endl;
}

}  // namespace Foam
