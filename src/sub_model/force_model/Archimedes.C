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

#include "./Archimedes.h"

namespace Foam {

cfdemDefineTypeName(Archimedes);

cfdemCreateNewFunctionAdder(forceModel, Archimedes);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
Archimedes::Archimedes(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      g_(cloud.globalF().g()) {
  createForceSubModels(subPropsDict_, kUnResolved);
  if (cloud_.modelType() == "A" || cloud_.modelType() == "Bfull") {
    CHECK(!forceSubModel_->treatForceBothCFDAndDEM())
        << ": model type A and Bfull request Archimedes model treatForceBothCFDAndDEM == false";
  } else if (cloud_.modelType() == "B") {
    CHECK(forceSubModel_->treatForceBothCFDAndDEM())
        << ": model type B request Archimedes model treatForceBothCFDAndDEM == true";
  }
}

Archimedes::~Archimedes() {}

void Archimedes::setForce() {
  base::MPI_Info("Setting Archimedes force...", true);
  Foam::vector buoyancy = Foam::vector::zero;
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkFineParticle(index) || cloud_.checkMiddleParticle(index)) {
      buoyancy = Foam::vector::zero;
      double raidus = cloud_.getRadius(index);
      // 这里使用 findCellID，不要使用 cellIDs，因为如果当前求解器不包含颗粒，则 cellIDs[index].mSize() == 0
      int findCellID = cloud_.findCellIDs()[index];
      if (findCellID >= 0) {
        double pV = 4.0 * M_PI * raidus * raidus * raidus / 3.0;
        buoyancy += -g_.value() * forceSubModel_->rhoField()[findCellID] * pV;
      }
      // write particle data to global array
      // index - particle index
      // buoyancy - total buoyancy
      forceSubModel_->partToArray(index, buoyancy, Foam::vector::zero, Foam::vector::zero, 0);
      if (forceSubModel_->verbose() && 0 == index) {
        Pout << "Archimedes buoyancy on particle " << index << ": [" << buoyancy[0] << ", " << buoyancy[1] << ", "
             << buoyancy[2] << "]" << endl;
      }
    }
  }
  base::MPI_Info("Setting Archimedes force - done", true);
}

}  // namespace Foam
