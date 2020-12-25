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
  Foam::vector particleVel = Foam::vector::zero;
  Foam::vector relativeVec = Foam::vector::zero;
  Foam::vector rotationVel = Foam::vector::zero;
  Foam::vector angularVel = Foam::vector::zero;
  Foam::vector updateVel = Foam::vector::zero;
  volVectorField globalUpdateVel = U;
  for (int index = 0; index < numberOfParticles(); ++index) {
    
  }  // End of loop all particles
}

}  // namespace Foam
