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
  The force model performs the calculation of forces (e.g. fluid-particle
  interaction forces) acting on each DEM particle. The ArchimedesIB model
  is a model that calculates the ArchimedesIB’s volumetric lift force
  stemming from density difference of fluid and particle. This model is
  especially suited for resolved CFD-DEM simulations where the particle
  is represented by immersed boundary method.

Syntax
  forceModels
  (
    ArchimedesIB
  );
  ArchimedesIBProps
  {
    gravityFieldName "g";
    voidfractionFieldName "voidfractionNext";
    treatForceExplicit true;
    verbose true;
  };

Restrictions
  Only for immersed boundary solvers.
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
      volumeFractionFieldName_(
          subPropsDict_.lookupOrDefault<Foam::word>("volumeFractionFieldName", "volumeFractionNext").c_str()),
      gravityFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("gravityFieldName", "g").c_str()),
      volumeFraction_(cloud.mesh().lookupObject<volScalarField>(volumeFractionFieldName_)),
#if defined(version21)
      g_(cloud.mesh().lookupObject<uniformDimensionedVectorField>(gravityFieldName_))
#elif defined(version16ext) || defined(version15)
      g_(dimensionedVector(
             cloud.mesh().lookupObject<IOdictionary>("environmentalProperties").lookup(environmentalProperties))
             .value())
#endif
{
  createForceSubModels(subPropsDict_, kResolved);
}

ArchimedesIB::~ArchimedesIB() {}

void ArchimedesIB::setForce() {
  Info << "Setting ArchimedesIB force..." << endl;
  Foam::vector buoyancy = Foam::vector::zero;
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // init
    buoyancy = Foam::vector::zero;
    // loop all mesh of current particle
    for (int subCell = 0; subCell < cloud_.particleOverMeshNumber()[index]; ++subCell) {
      label cellI = cloud_.cellIDs()[index][subCell];
      if (cellI > -1) {  // cell found
        buoyancy +=
            -g_.value() * forceSubModel_->rhoField()[cellI] * cloud_.mesh().V()[cellI] * (1.0 - volumeFraction_[cellI]);
      }
    }
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
