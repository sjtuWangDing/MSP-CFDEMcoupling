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

#include "./implicit_couple.h"
#include "sub_model/data_exchange_model/data_exchange_model.h"
#include "sub_model/force_model/global_force.h"

namespace Foam {

cfdemDefineTypeName(implicitCouple);

cfdemCreateNewFunctionAdder(momCoupleModel, implicitCouple);

//! \brief Constructor
implicitCouple::implicitCouple(cfdemCloud& cloud)
    : momCoupleModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      velFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("velFieldName", "U").c_str()),
      granVelFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("granVelFieldName", "Us").c_str()),
      voidFractionFieldName_(
          subPropsDict_.lookupOrDefault<Foam::word>("voidFractionFieldName", "voidFraction").c_str()),
      U_(cloud.mesh().lookupObject<volVectorField>(velFieldName_)),
      Us_(cloud.mesh().lookupObject<volVectorField>(granVelFieldName_)),
      alpha_(cloud.mesh().lookupObject<volScalarField>(voidFractionFieldName_)),
      KslLimit_(1e10),
      KslPrev_(IOobject("KslPrev", cloud.mesh().time().timeName(), cloud.mesh(),
                        IOobject::READ_IF_PRESENT,  // MUST_READ,
                        IOobject::NO_WRITE),
               cloud.mesh().lookupObject<volScalarField>("Ksl")),
      KslNext_(IOobject("KslNext", cloud.mesh().time().timeName(), cloud.mesh(),
                        IOobject::READ_IF_PRESENT,  // MUST_READ,
                        IOobject::NO_WRITE),
               cloud.mesh().lookupObject<volScalarField>("Ksl")) {
  if (subPropsDict_.found("KslLimit")) {
    KslLimit_ = readScalar(subPropsDict_.lookup("KslLimit"));
    Info << "implicit momentum exchange field is limited to : " << KslLimit_ << endl;
  }
  if (subPropsDict_.found("maxAlpha")) {
    maxAlpha_ = readScalar(subPropsDict_.lookup("maxAlpha"));
    Info << "implicit momentum exchange field calculate if fluid volume fraction smaller than : " << maxAlpha_ << endl;
  }
}

//! \brief Destructor
implicitCouple::~implicitCouple() {}

//! \brief reset mom field
void implicitCouple::resetMomSourceField() {
  KslPrev_ == KslNext_;
  KslNext_ == dimensionedScalar("zero", KslNext_.dimensions(), 0.0);
}

//! \brief get implicit momentum source field
tmp<volScalarField> implicitCouple::impMomSource() {
  // update KslNext in first sub-time step
  if (cloud_.dataExchangeM().checkValidCouplingStep()) {
    forAll(KslNext_, cellI) {
      double magUr = mag(U_[cellI] - Us_[cellI]);
      double magImplForce = mag(cloud_.globalF().impParticleForce()[cellI]);
      if (magUr > Foam::SMALL && magImplForce > Foam::SMALL && alpha_[cellI] < maxAlpha_) {
        KslNext_[cellI] = magImplForce / magUr / cloud_.mesh().V()[cellI];
      } else {
        KslNext_[cellI] = 0;
      }
      if (KslNext_[cellI] > KslLimit_) {
        KslNext_[cellI] = KslLimit_;
      }
    }
  }
  double tsf = cloud_.dataExchangeM().timeStepFraction();
  return tmp<volScalarField>(new volScalarField("KslImplicitCouple", tsf * KslPrev_ + (1 - tsf) * KslNext_));
}

}  // namespace Foam
