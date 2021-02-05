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
#include "cloud/cfdem_cloud_semi.h"
#include "cfdem_tools/cfdem_tools.h"

#if defined(version30)
  #include "turbulentTransportModel.H"
  #include "pisoControl.H"
#else
  #include "turbulenceModel.H"
#endif

#if defined(versionv1606plus) || defined(version40)
  #include "fvOptions.H"
#else
  #include "fvIOoptionList.H"
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
  #include "createFvOptions.H"
  #include "initContinuityErrs.H"
  #include "readGravitationalAcceleration.H"
  // create cfdemCloud
  Foam::cfdemCloudSemi particleCloud(mesh);
  cfdemTools::checkModelType(particleCloud);
  std::string modelType = particleCloud.modelType();

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
    particleCloud.evolve(U, voidFraction, Us, Ksl);

    surfaceScalarField voidFractionf = fvc::interpolate(voidFraction);
    phi = voidFractionf * phiByVoidFraction;

    // #include "forceCheckIm.H"
    Info << "\nSolver level total Eulerian momentum exchange:"<< endl;
    volVectorField fImp(Ksl * (Us - U));
    particleCloud.scaleWithVcell(fImp);
    dimensionedVector fImpTotal = gSum(fImp);
    Info << "  TotalForceImp:  " << fImpTotal.value() << endl;
    Info << "  Warning, these values are based on latest Ksl and Us but prev. iteration U!\n" << endl;

    if (particleCloud.solveFlow()) {
      // Momentum predictor
      fvVectorMatrix UEqn
      (
          fvm::ddt(voidFraction, U) - fvm::Sp(fvc::ddt(voidFraction), U)
        + fvm::div(phi, U) - fvm::Sp(fvc::div(phi), U)
        + particleCloud.divVoidFractionTau(U, voidFraction)
        ==
        - fvm::Sp(Ksl / rho, U)
        + fvOptions(U)
      );
      UEqn.relax();
      fvOptions.constrain(UEqn);

    #if defined(version30)
      if (piso.momentumPredictor())
    #else
      if (momentumPredictor)
    #endif
      {
        if (modelType == "B" || modelType == "Bfull") {
          solve(UEqn == - fvc::grad(p) + Ksl / rho * Us);
        } else {
          solve(UEqn == - voidFraction * fvc::grad(p) + Ksl / rho * Us);
        }
        fvOptions.correct(U);
      }
    // PISO loop
    #if defined(version30)
      while (piso.correct())
    #else
      for (int corr = 0; corr < nCorr; corr++)
    #endif
      {
        volScalarField rUA = 1.0 / UEqn.A();
        surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));
        volScalarField rUAvoidFraction("(voidfraction2|A(U))",rUA * voidFraction);
        surfaceScalarField rUAfvoidFraction("(voidfraction2|A(U)F)", fvc::interpolate(rUAvoidFraction));
        U = rUA * UEqn.H();
      #ifdef version23
        phi = (fvc::interpolate(U) & mesh.Sf()) + rUAfvoidFraction * fvc::ddtCorr(U, phiByVoidFraction);
      #else
        phi = (fvc::interpolate(U) & mesh.Sf()) + fvc::ddtPhiCorr(rUAvoidFraction, U, phiByVoidFraction);
      #endif
        surfaceScalarField phiS(fvc::interpolate(Us) & mesh.Sf());
        phi += rUAf * (fvc::interpolate(Ksl / rho) * phiS);

        if (modelType=="A") {
          rUAvoidFraction = volScalarField("(voidfraction2|A(U))",rUA * voidFraction * voidFraction);
        }

        #include "./fixedFluxPressureHandling.h"

        // Non-orthogonal pressure corrector loop
      #if defined(version30)
        while (piso.correctNonOrthogonal())
      #else
        for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
      #endif
        {
          // Pressure corrector
          fvScalarMatrix pEqn
          (
            fvm::laplacian(rUAvoidFraction, p) == fvc::div(voidFractionf * phi) + particleCloud.ddtVoidFraction()
          );
          pEqn.setReference(pRefCell, pRefValue);

        #if defined(version30)
          pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
          if (piso.finalNonOrthogonalIter()) {
            phiByVoidFraction = phi - pEqn.flux() / voidFractionf;
          }
        #else
          if (corr == nCorr-1 && nonOrth == nNonOrthCorr) {
          #if defined(versionExt32)
            pEqn.solve(mesh.solutionDict().solver("pFinal"));
          #else
            pEqn.solve(mesh.solver("pFinal"));
          #endif
          } else {
            pEqn.solve();
          }
          if (nonOrth == nNonOrthCorr) {
            phiByVoidFraction = phi - pEqn.flux() / voidFractionf;
          }
        #endif
        }  // end non-orthogonal corrector loop
        phi = voidFractionf * phiByVoidFraction;
        #include "./continuityErrorPhiPU.h"

        if (modelType == "B" || modelType == "Bfull") {
          U -= rUA * fvc::grad(p) - Ksl / rho * Us * rUA;
        } else {
          U -= voidFraction * rUA * fvc::grad(p) - Ksl / rho * Us * rUA;
        }
        U.correctBoundaryConditions();
        fvOptions.correct(U);
      }  // end piso loop
      laminarTransport.correct();
      turbulence->correct();
    }  // solveFlow

    runTime.write();
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << ", ClockTime = " << runTime.elapsedClockTime() << " s" << endl;
  }  // end of runtime loop
  Info << "cfdemSolverPiso - done\n" << endl;
}
