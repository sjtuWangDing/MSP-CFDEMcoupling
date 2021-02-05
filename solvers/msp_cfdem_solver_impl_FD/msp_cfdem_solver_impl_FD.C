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
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "pisoControl.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "dynamicFvMesh.H"
#include "cloud/cfdem_cloud_IB_opti.h"

int main(int argc, char* argv[]) {
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createDynamicFvMesh.H"
  #include "createControl.H"
  #include "createFields.H"
  #include "createFvOptions.H"
  #include "initContinuityErrs.H"

  turbulence->validate();

  // create cfdemCloud
  Foam::cfdemCloudIBOpti particleCloud(mesh);

  Info << "\nStarting time loop\n" << endl;

  while(runTime.loop()) {

    Info << "Time = " << runTime.timeName() << nl << endl;

    #include "CourantNo.H"

    // particle cloud evolve and update mesh
    particleCloud.evolve(volumeFraction, interface);

    if (particleCloud.solveFlow()) {
      // Momentum predictor
      fvVectorMatrix UEqn
      (
          fvm::ddt(U)
        + fvm::div(phi, U)
        + turbulence->divDevReff(U)
        ==
          fvOptions(U)
      );
      UEqn.relax();
      fvOptions.constrain(UEqn);

      if (piso.momentumPredictor()) {
        solve(UEqn == -fvc::grad(p));
        fvOptions.correct(U);
      }

      // PISO loop
      while (piso.correct()) {
        // rAU == 1 / Ad，是标量场，且已知，对应于线性方程组主对角线系数
        volScalarField rAU(1.0 / UEqn.A());

        // HbyA = Ah / Ad，是矢量场，且已知，因为 Ah = As - An * U，这里 U 是预测速度(未修正)
        // volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));
        volVectorField HbyA("HbyA", U);
        HbyA = rAU * UEqn.H();

        // 定义 HbyA 的通量，是面标量场，且已知
        // HbyA = fvc::interpolate(HbyA) & mesh.Sf() + correction(修正项)
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            fvc::flux(HbyA)  // fvc::interpolate(HbyA) & mesh.Sf()
          + fvc::interpolate(rAU) * fvc::ddtCorr(U, phi)  // correction
        );

        adjustPhi(phiHbyA, U, p);

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p, U, phiHbyA, rAU);

        // Non-orthogonal pressure corrector loop
        while (piso.correctNonOrthogonal()) {
          // p 方程
          fvScalarMatrix pEqn
          (
            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
          );
          pEqn.setReference(pRefCell, pRefValue);
          // pressure corrector
          pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
          if (piso.finalNonOrthogonalIter()) {
            // 更新通量场 phi
            // U = HbyA - grad(p) / Ad  ===>  phi = (HbyA - grad(p)) & mesh.Sf()
            // pEqn.flux() 返回的就是 grad(p) & mesh.Sf()
            phi = phiHbyA - pEqn.flux();
          }
        }
        #include "continuityErrs.H"
        // 通过新的压力场 p 更新速度场 U
        U = HbyA - rAU * fvc::grad(p);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
      } // piso loop
    } // solve flow

    // 针对非牛顿流体的修正
    laminarTransport.correct();
    // 如果指定了湍流模型，则解额外的湍流方程
    turbulence->correct();

    // 通过颗粒速度以及phiIB修正速度与压力
    volScalarField volumeFractionNext = mesh.lookupObject<volScalarField>("volumeFractionNext");
    particleCloud.calcVelocityCorrection(p, U, phiIB, volumeFractionNext);
    fvOptions.correct(U);

    runTime.write();
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << " ClockTime = " << runTime.elapsedClockTime() << " s"
         << endl;
  } // end of runtime loop

  Info << "mspCfdemSolverImplFD - done\n" << endl;
  return 0;
}
