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

#include "./Archimedes_IB_test.h"

namespace Foam {

cfdemDefineTypeName(ArchimedesIBTest);

cfdemCreateNewFunctionAdder(forceModel, ArchimedesIBTest);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
ArchimedesIBTest::ArchimedesIBTest(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      g_(cloud.globalF().g()) {
  createForceSubModels(subPropsDict_, kResolved);
}

ArchimedesIBTest::~ArchimedesIBTest() {}

void ArchimedesIBTest::setForce() {
  Info << "Setting ArchimedesIBTest force..." << endl;
  Info << "Warning: ArchimedesIBTest::setForce(): "
       << "directly use particle's buoyancy force with rho of fluid == 1 and number of processors == 4" << endl;
  base::MPI_Barrier(0.5);
  Foam::vector buoyancy = Foam::vector::zero;
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    buoyancy = Foam::vector::zero;
    scalar radius = cloud_.getRadius(index);
    scalar volume = (4.0 / 3.0) * M_PI * radius * radius * radius / 4.0;
    scalar rho = 1;
    buoyancy = -g_.value() * rho * volume;
    Pout << "force: " << buoyancy[0] << ", " << buoyancy[1] << ", " << buoyancy[2] << endl;
    base::MPI_Barrier(0.5);
    // write particle data to global array
    // index - particle index
    // buoyancy - total buoyancy
    forceSubModel_->partToArray(index, buoyancy, Foam::vector::zero, Foam::vector::zero, 0);
    if (forceSubModel_->verbose()) {
      Info << "Archimedes buoyancy on particle " << index << ": [" << buoyancy[0] << ", " << buoyancy[1] << ", "
           << buoyancy[2] << "]" << endl;
    }
  }
}

}  // namespace Foam
