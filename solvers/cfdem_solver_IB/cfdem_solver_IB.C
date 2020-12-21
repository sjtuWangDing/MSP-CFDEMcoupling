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
  cfdemSolverIB is a coupled CFD-DEM solver using CFDEMcoupling, an open
  source parallel coupled CFD-DEM framework, for calculating the dynamics
  between immersed bodies and the surrounding fluid. Being an implementation
  of an immersed boundary method it allows tackling problems where the body
  diameter exceeds the maximal size of a fluid cell. Using the toolbox of
  OpenFOAM®(*) the governing equations of the fluid are computed and the
  corrections of velocity and pressure field with respect to the body-
  movement information, gained from LIGGGHTS, are incorporated.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"

#include "cloud/of_version.h"
#include "cloud/cfdem_cloud_IB.h"

#if defined(version30)
  #include "turbulentTransportModel.H"
  #include "pisoControl.H"
#else
  #include "turbulenceModel.H"
#endif

int main(int argc, char* argv[]) {

  #include "setRootCase.H"

  #include "createTime.H"

  #include "createDynamicFvMesh.H"

  #if defined(version30)
    pisoControl piso(mesh);
    #include "createTimeControls.H"
  #endif

  #include "create_fields.h"

  #include "initContinuityErrs.H"

  #if defined(version22)
    #include "createFvOptions.H"
  #endif

  // create cfdemCloud
  Foam::cfdemCloudIB particleCloud(mesh);

  Info << "\nStarting time loop\n" << endl;

  while(runTime.loop()) {

    Info << "Time = " << runTime.timeName() << endl << endl;

    // set interface
    interface = mag(mesh.lookupObject<volScalarField>("volumeFractionNext"));
    particleCloud.setMeshHasUpdated(mesh.update());

    #if defined(version30)
      #include "readTimeControls.H"
      #include "CourantNo.H"
      #include "setDeltaT.H"
    #else
      #include "readPISOControls.H"
      #include "CourantNo.H"
    #endif

    // particle evolve
    particleCloud.evolve(volumeFraction, interface);

  } // end of runtime loop

  Info << "cfdemCloudIB - done\n" << endl;

  return 0;
}
