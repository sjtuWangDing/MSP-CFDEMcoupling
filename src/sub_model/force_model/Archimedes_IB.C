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

#include "./Archimedes_IB.h"

namespace Foam {

cfdemDefineTypeName(ArchimedesIB);

cfdemCreateNewFunctionAdder(forceModel, ArchimedesIB);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
ArchimedesIB::ArchimedesIB(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      volumeFraction_(cloud_.globalF().volumeFraction()),
      g_(cloud_.globalF().g()) {
  createForceSubModels(subPropsDict_, kResolved);
}

ArchimedesIB::~ArchimedesIB() {}

void ArchimedesIB::setForce() {
  Info << "Setting ArchimedesIB force..." << endl;
  Foam::vector buoyancy = Foam::vector::zero;
  std::once_flag onceOp;
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkCoarseParticle(index)) {
      // init
      buoyancy = Foam::vector::zero;
      // loop all mesh of current particle
      for (int subCell = 0; subCell < cloud_.particleOverMeshNumber()[index]; ++subCell) {
        int cellID = cloud_.cellIDs()[index][subCell];
        if (cellID >= 0) {  // cell found
          buoyancy += -g_.value() * forceSubModel_->rhoField()[cellID] * cloud_.mesh().V()[cellID] *
                      (1.0 - volumeFraction_[cellID]);
        }
      }
      // write particle data to global array
      // index - particle index
      // buoyancy - total buoyancy
      forceSubModel_->partToArray(index, buoyancy, Foam::vector::zero, Foam::vector::zero, 0);
      std::call_once(onceOp, [this, &buoyancy, &index]() {
        if (forceSubModel_->verbose() && mag(buoyancy) > Foam::SMALL) {
          Pout << "index = " << index << endl;
          Pout << "buoyancy (part) = " << buoyancy << endl;
        }
      });
    }
  }
  Info << "Setting ArchimedesIB force - done" << endl;
}

}  // namespace Foam
