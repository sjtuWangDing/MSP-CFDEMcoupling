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
  cfdemCloudIBLmpf derived from Foam::cfdemCloudIBOpti

Class
  Foam::cfdemCloudIBLmpf
\*---------------------------------------------------------------------------*/

#include "cloud/cfdem_cloud_IB_lmpf.h"

namespace Foam {

//! \brief Constructed from mesh
cfdemCloudIBLmpf::cfdemCloudIBLmpf(const fvMesh& mesh) : cfdemCloudIBOpti(mesh) {}

//! \brief Destructor
cfdemCloudIBLmpf::~cfdemCloudIBLmpf() {}

void cfdemCloudIBLmpf::calcLmpf(const volVectorField& U, const volScalarField& rho,
                                const volScalarField& volumeFraction, volVectorField& lmpf) const {
  // set particle velocity
  Foam::vector parPos = Foam::vector::zero;
  Foam::vector lVel = Foam::vector::zero;
  Foam::vector angVel = Foam::vector::zero;
  Foam::vector rVec = Foam::vector::zero;
  Foam::vector parVel = Foam::vector::zero;
  // reset lmpf field
  lmpf == dimensionedVector("zero", lmpf.dimensions(), vector::zero);
  // 计算更新速度
  volVectorField globalUpdateVel = U;
  for (int index = 0; index < numberOfParticles(); ++index) {
    parPos = getPosition(index);         // 颗粒中心
    lVel = getVelocity(index);           // 颗粒线速度
    angVel = getAngularVelocity(index);  // 颗粒角速度
    for (int subCell = 0; subCell < particleOverMeshNumber()[index]; ++subCell) {
      int cellID = cellIDs()[index][subCell];
      if (cellID > -1) {
        for (int i = 0; i < 3; ++i) {
          rVec[i] = U.mesh().C()[cellID][i] - parPos[i];
        }
        // 计算颗粒速度
        parVel = lVel + (angVel ^ rVec);
        double vf = volumeFractions()[index][subCell];
        globalUpdateVel[cellID] = (1.0 - vf) * parVel + vf * U[cellID];
      }
    }  // End of loop all meshes
  }    // End of loop all particles
  Foam::vector sumLmpf = Foam::vector::zero;
  Foam::vector sumVelDiff = Foam::vector::zero;
  for (int index = 0; index < numberOfParticles(); ++index) {
    for (int subCell = 0; subCell < particleOverMeshNumber()[index]; ++subCell) {
      int cellID = cellIDs()[index][subCell];
      if (cellID > -1) {
        Foam::vector updateVel = Foam::vector::zero;
        Foam::vector neighUpdateVel = Foam::vector::zero;
        // 获取相邻网格的所有 neighbour cell
        const labelList& neighList = U.mesh().cellCells()[cellID];
        int neighbourNum = 0;
        forAll(neighList, i) {
          int neighbourCellID = neighList[i];
          if (neighbourCellID < 0) {
            continue;
          }
          neighbourNum += 1;
          neighUpdateVel += globalUpdateVel[neighbourCellID];
        }  // End of neighbour loop
        if (neighbourNum > 0) {
          updateVel = 0.5 * globalUpdateVel[cellID] + 1.0 * neighUpdateVel / (2.0 * neighbourNum);
        } else {
          updateVel = globalUpdateVel[cellID];
        }
        // 计算 lmpf
        lmpf[cellID] = (updateVel - U[cellID]) / (U.mesh().time().deltaT().value());
        sumLmpf += lmpf[cellID] * rho[cellID] * U.mesh().V()[cellID];
        sumVelDiff += updateVel - U[cellID];
      }
    }  // End of loop all meshes
  }    // End of loop all particles
  if (mag(sumVelDiff) > Foam::SMALL) {
    Pout << "sumVelDiff: " << sumVelDiff << endl;
  }
  base::MPI_Barrier(0.1);
  if (mag(sumLmpf) > Foam::SMALL) {
    Pout << "sumLmpf: " << sumLmpf << endl;
  }
  base::MPI_Barrier(0.1);
}

}  // namespace Foam
