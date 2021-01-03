/*---------------------------------------------------------------------------*\
  CFDEMcoupling - Open Source CFD-DEM coupling

  CFDEMcoupling is part of the CFDEMproject
  www.cfdem.com
                              Christoph Goniva, christoph.goniva@cfdem.com
                              Copyright (C) 1991-2009 OpenCFD Ltd.
                              Copyright (C) 2009-2012 JKU, Linz
                              Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
  This file is part of CFDEMcoupling.

  CFDEMcoupling is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Description
  cfdemSolverPiso is a coupled CFD-DEM solver using CFDEMcoupling,
  an open source parallel coupled CFD-DEM framework. Based on pisoFoam,
  a finite volume based solver for turbulent Navier-Stokes equations applying
  the PISO algorithm, cfdemSolverPiso has additional functionality for
  a coupling to the DEM code LIGGGHTS. The volume averaged Navier-Stokes
  Equations are solved accounting for momentum exchange and volume
  displacement of discrete particles whose trajectories are calculated
  in the DEM code LIGGGHTS.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"

#include "cloud/of_version.h"
#include "cloud/cfdem_cloud.h"

#if defined(version30)
  #include "turbulentTransportModel.H"
  #include "pisoControl.H"
#else
  #include "turbulenceModel.H"
#endif

int main(int argc, char *argv[]) {
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  #if defined(version30)
    pisoControl piso(mesh);
    #include "createTimeControls.H"
  #endif
  #include "create_fields.h"
  #include "initContinuityErrs.H"
  #include "readGravitationalAcceleration.H"
  // create cfdemCloud
  Foam::cfdemCloud particleCloud(mesh);
  // run loop
  Info << "\nStarting time loop\n" << endl;
  while (runTime.loop()) {
    Info << "\nTime = " << runTime.timeName() << endl << endl;
    #if defined(version30)
      #include "readTimeControls.H"
      #include "CourantNo.H"
      #include "setDeltaT.H"
    #else
      #include "readPISOControls.H"
      #include "CourantNo.H"
    #endif
    particleCloud.evolve(voidFraction, Us, U);
  }
}
