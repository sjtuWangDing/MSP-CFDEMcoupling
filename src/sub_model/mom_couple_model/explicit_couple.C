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

#include "./explicit_couple.h"
#include "sub_model/data_exchange_model/data_exchange_model.h"
#include "sub_model/force_model/global_force.h"

namespace Foam {

cfdemDefineTypeName(explicitCouple);

cfdemCreateNewFunctionAdder(momCoupleModel, explicitCouple);

//! \brief Constructor
explicitCouple::explicitCouple(cfdemCloud& cloud)
    : momCoupleModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      expForceLimit_(1.0e10, 1.0e10, 1.0e10),
      expForcePrev_(IOobject("expForcePrev", cloud.mesh().time().timeName(), cloud.mesh(),
                             IOobject::READ_IF_PRESENT,  // MUST_READ,
                             IOobject::NO_WRITE),
                    cloud.mesh().lookupObject<volVectorField>("expForce")),
      expForceNext_(IOobject("expForceNext", cloud.mesh().time().timeName(), cloud.mesh(),
                             IOobject::READ_IF_PRESENT,  // MUST_READ,
                             IOobject::NO_WRITE),
                    cloud.mesh().lookupObject<volVectorField>("expForce")) {
  if (subPropsDict_.found("expForceLimit")) {
    expForceLimit_ = Foam::vector(subPropsDict_.lookup("expForceLimit"));
    Info << "explicit momentum exchange field is limited to : " << expForceLimit_ << endl;
  }
}

//! \brief Destructor
explicitCouple::~explicitCouple() {}

//! \brief reset mom field
void explicitCouple::resetMomSourceField() {
  expForcePrev_ == expForceNext_;
  expForceNext_ == dimensionedVector("zero", expForceNext_.dimensions(), Foam::vector::zero);
}

//! \brief get explicit momentum source field
tmp<volVectorField> explicitCouple::expMomSource() {
  // update expForceNext_ in first sub-time step
  if (cloud_.dataExchangeM().checkValidCouplingStep()) {
    forAll(expForceNext_, cellI) {
      expForceNext_[cellI] = cloud_.globalF().expParticleForce()[cellI] / cloud_.mesh().V()[cellI];
#pragma unroll
      for (int i = 0; i < 3; ++i) {
        if (expForceNext_[cellI][i] > expForceLimit_[i]) {
          expForceNext_[cellI][i] = expForceLimit_[i];
        }
      }
    }
  }
  double tsf = cloud_.dataExchangeM().timeStepFraction();
  return tmp<volVectorField>(new volVectorField("expForceCouple", (1 - tsf) * expForcePrev_ + tsf * expForceNext_));
}

}  // namespace Foam
