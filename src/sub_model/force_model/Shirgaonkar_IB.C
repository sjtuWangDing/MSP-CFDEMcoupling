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
      fictitiousForce_(cloud.globalF().fictitiousForce()),
      volumeFraction_(cloud.globalF().volumeFraction()),
      useTorque_(subPropsDict_.lookupOrDefault<bool>("useTorque", false)) {
  createForceSubModels(subPropsDict_, kResolved);
}

ShirgaonkarIB::~ShirgaonkarIB() {}

void ShirgaonkarIB::setForce() {
  base::MPI_Info("Setting ShirgaonkarIB force...", true);
  volVectorField IBDrag = forceSubModel_->IBDrag(U_, p_);
  Foam::vector particleCenterPos = Foam::vector::zero;
  Foam::vector cellPos = Foam::vector::zero;
  Foam::vector drag = Foam::vector::zero;
  Foam::vector torque = Foam::vector::zero;
  Foam::vector particleU = Foam::vector::zero;
  double radius = 0.0;
  std::once_flag onceOp;
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkCoarseParticle(index)) {
      // init
      drag = Foam::vector::zero;
      torque = Foam::vector::zero;
      // get index's particle center position
      particleCenterPos = cloud_.getPosition(index);
      particleU = cloud_.getVelocity(index);
      radius = cloud_.getRadius(index);
      // loop all mesh of current particle
      for (int subCell = 0; subCell < cloud_.particleOverMeshNumber()[index]; ++subCell) {
        int cellID = cloud_.cellIDs()[index][subCell];
        if (cellID >= 0) {  // cell Found
          cellPos = cloud_.mesh().C()[cellID];
#if 1
          // 计算颗粒雷诺数
          double rhoCell = cloud_.globalF().rhoField()[cellID];
          double nuf = forceSubModel_->nuField()[cellID];
          double vCell = cloud_.mesh().V()[cellID];
          double Re = mag(particleU) * 2 * radius / (nuf + Foam::SMALL);
          double vf = volumeFraction_[cellID];
          drag += IBDrag[cellID] * IBDrag.mesh().V()[cellID] * (1 - vf);
          torque += (cellPos - particleCenterPos) ^ IBDrag[cellID] * (1 - vf) * IBDrag.mesh().V()[cellID];
          if (Re < 50) {
            drag -= fictitiousForce_[cellID] * rhoCell * vCell;
            torque -= (cellPos - particleCenterPos) ^ fictitiousForce_[cellID] * rhoCell * vCell;
          } else {
            drag -= fictitiousForce_[cellID] * rhoCell * vCell * (1 - vf);
            torque -= (cellPos - particleCenterPos) ^ fictitiousForce_[cellID] * rhoCell * vCell * (1 - vf);
          }
#else
          drag += IBDrag[cellID] * IBDrag.mesh().V()[cellID];
          torque += (cellPos - particleCenterPos) ^ IBDrag[cellID] * IBDrag.mesh().V()[cellID];
#endif
        }
      }
      // write particle data to global array
      forceSubModel_->partToArray(index, drag, Foam::vector::zero, Foam::vector::zero, 0);

      std::call_once(onceOp, [this, &drag, &index]() {
        if (forceSubModel_->verbose() && mag(drag) > Foam::SMALL) {
          Pout << "index = " << index << endl;
          Pout << "drag (part) = " << drag << endl;
        }
      });

      if (useTorque_) {
        forceSubModel_->addTorque(index, torque);
      }
    }
  }
  base::MPI_Info("Setting ShirgaonkarIB force - done", true);
}

}  // namespace Foam
