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

    if (particleCloud.solveFlow()) {
      // Momentum predictor
      fvVectorMatrix UEqn(
          fvm::ddt(U)
        + fvm::div(phi, U)
        + turbulence->divDevReff(U)
        #if defined(version22)
        ==
          fvOptions(U)
        #endif
      );

      UEqn.relax();
      #if defined(version22)
      fvOptions.constrain(UEqn);
      #endif

      #if defined(version30)
      if (piso.momentumPredictor())
      #else
      if (momentumPredictor)
      #endif
      {
        solve(UEqn == -fvc::grad(p));
      }

      // PISO loop
      #if defined(version30)
      while (piso.correct())
      #else
      for (int corr = 0; corr < nCorr; corr++)
      #endif
      {
        volScalarField rUA = 1.0 / UEqn.A();
        surfaceScalarField rUAf(fvc::interpolate(rUA));
        U = rUA * UEqn.H();
        #ifdef version23
        phi = (fvc::interpolate(U) & mesh.Sf()) + rUAf * fvc::ddtCorr(U, phi);
        #else
        phi = (fvc::interpolate(U) & mesh.Sf()) + fvc::ddtPhiCorr(rUA, U, phi);
        #endif
        adjustPhi(phi, U, p);

        #if defined(version22)
        fvOptions.relativeFlux(phi);
        #endif

        // Non-orthogonal pressure corrector loop
        #if defined(version30)
        while (piso.correctNonOrthogonal())
        #else
        for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
        #endif
        {
          // Pressure corrector
          fvScalarMatrix pEqn(
            fvm::laplacian(rUA, p) == fvc::div(phi) + particleCloud.ddtVoidFraction()
          );
          pEqn.setReference(pRefCell, pRefValue);
          #if defined(version30)
          pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
          if (piso.finalNonOrthogonalIter()) {
            phi -= pEqn.flux();
          }
          #else
          if (corr == nCorr - 1 && nonOrth == nNonOrthCorr) {
            #if defined(versionExt32)
            pEqn.solve(mesh.solutionDict().solver("pFinal"));
            #else
            pEqn.solve(mesh.solver("pFinal"));
            #endif
          } else {
            pEqn.solve();
          }
          if (nonOrth == nNonOrthCorr) {
            phi -= pEqn.flux();
          }
          #endif
        } // pressure corrector
        #include "continuityErrs.H"
        U -= rUA * fvc::grad(p);
        U.correctBoundaryConditions();
      } // piso loop
    } // solve flow

    laminarTransport.correct();
    turbulence->correct();
    volScalarField volumeFractionNext = mesh.lookupObject<volScalarField>("volumeFractionNext");
    particleCloud.calcVelocityCorrection(p, U, phiIB, volumeFractionNext);
    #if defined(version22)
    fvOptions.correct(U);
    #endif
    runTime.write();
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << " ClockTime = " << runTime.elapsedClockTime() << " s"
         << endl;
  } // end of runtime loop

  Info << "cfdemCloudIB - done\n" << endl;

  return 0;
}
