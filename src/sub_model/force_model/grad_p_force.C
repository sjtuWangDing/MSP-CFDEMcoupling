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

#include "./grad_p_force.h"

namespace Foam {

cfdemDefineTypeName(gradPForce);

cfdemCreateNewFunctionAdder(forceModel, gradPForce);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
gradPForce::gradPForce(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      pressureFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("pressureFieldName", "p").c_str()),
      p_(cloud.mesh().lookupObject<volScalarField>(pressureFieldName_)) {
  createForceSubModels(subPropsDict_, kResolved);
  if (cloud_.modelType() == "B") {
    FatalError << "using gradPForce with model type B is not valid\n" << abort(FatalError);
  }
  if (!forceSubModel_->treatForceExplicitInMomEquation()) {
    FatalError << "gradPForce model need treatForceExplicitInMomEquation = true\n" << abort(FatalError);
  }
  if (cloud_.modelType() == "A" && forceSubModel_->treatForceBothCFDAndDEM()) {
    FatalError << "gradPForce with model type A requires treatForceBothCFDAndDEM = false\n" << abort(FatalError);
  }
  if (cloud_.modelType() == "Bfull" && !forceSubModel_->treatForceBothCFDAndDEM()) {
    FatalError << "gradPForce with model type Bfull requires treatForceBothCFDAndDEM = true\n" << abort(FatalError);
  }
}

gradPForce::~gradPForce() {}

void gradPForce::setForce() {
  Info << "Setting grad p force..." << endl;
  volVectorField gradP_ = fvc::grad(p_);
  gradPInterpolator_.reset(
      interpolation<vector>::New(subPropsDict_.lookupOrDefault("gradPInterpolationType", word("cellPointFace")), gradP_)
          .ptr());
  double radius = 0.0;
  Foam::vector force = Foam::vector::zero;
  Foam::vector gradP = Foam::vector::zero;
  Foam::vector particlePos = Foam::vector::zero;
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    radius = cloud_.getRadius(index);
    particlePos = cloud_.getPosition(index);
    force = Foam::vector::zero;
    int findCellID = cloud_.findCellIDs()[index];
    if (findCellID >= 0) {
      if (forceSubModel_->interpolation()) {
        gradP = gradPInterpolator_().interpolate(particlePos, findCellID);
      } else {
        gradP = gradP_[findCellID];
      }
      force = -gradP * forceSubModel_->rhoField()[findCellID] * cloud_.voidFractionM().pV(radius);
    }
    if (forceSubModel_->verbose()) {
      Pout << "grad p force: " << force << endl;
      Pout << "Archimedes force: " << -9.81 * forceSubModel_->rhoField()[findCellID] * cloud_.voidFractionM().pV(radius)
           << endl;
    }
    // write particle data to global array
    forceSubModel_->partToArray(index, force, Foam::vector::zero, Foam::vector::zero, 0.0);
  }
  Info << "Setting grad p force..." << endl;
}

}  // namespace Foam
