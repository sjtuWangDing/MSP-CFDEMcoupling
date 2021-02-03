/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pisoFoam

Description
    Transient solver for incompressible, turbulent flow, using the PISO
    algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "pisoControl.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H" // 湍流模型头文件，在 create_fields.h 文件中定义了湍流模型的对象

int main(int argc, char* argv[]) {
  #include "postProcess.H"
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  #include "createControl.H"
  #include "createFields.H"
  #include "createFvOptions.H"
  #include "initContinuityErrs.H"

  turbulence->validate();

  Info << "\nStarting time loop\n" << endl;

  while (runTime.loop()) {
    Info << "Time = " << runTime.timeName() << nl << endl;
    #include "CourantNo.H"

    // Pressure-velocity PISO corrector
    // Solve the Momentum equation
    MRF.correctBoundaryVelocity(U);
    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + MRF.DDt(U)  // relates to rotating framework(旋转坐标系)
      + turbulence->divDevReff(U)  // 考虑了湍流模型的扩散项，所以这其中不仅包含了 nu，还有比如 Boussinesq 假设引入的扩散系数 nut
      ==
        fvOptions(U)  // 可以添加自定义的源项
    );
    UEqn.relax();
    fvOptions.constrain(UEqn);

    if (piso.momentumPredictor()) {
      solve(UEqn == -fvc::grad(p));
      fvOptions.correct(U);
    }

    // PISO loop
    while (piso.correct()) {
      volScalarField rAU(1.0 / UEqn.A());
      volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));
      surfaceScalarField phiHbyA
      (
        "phiHbyA",
        fvc::flux(HbyA) + fvc::interpolate(rAU) * fvc::ddtCorr(U, phi)
      );

      MRF.makeRelative(phiHbyA);

      adjustPhi(phiHbyA, U, p);

      // Update the pressure BCs to ensure flux consistency
      constrainPressure(p, U, phiHbyA, rAU, MRF);

      // Non-orthogonal pressure corrector loop
      while (piso.correctNonOrthogonal()) {
          // Pressure corrector
          fvScalarMatrix pEqn
          (
            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
          );

          pEqn.setReference(pRefCell, pRefValue);

          pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

          if (piso.finalNonOrthogonalIter()) {
            phi = phiHbyA - pEqn.flux();
          }
      }

      #include "continuityErrs.H"

      U = HbyA - rAU * fvc::grad(p);
      U.correctBoundaryConditions();
      fvOptions.correct(U);
    }

    laminarTransport.correct(); // 针对非牛顿流体的修正
    turbulence->correct(); // 如果指定了湍流模型，比如 k-w，那么是需要求解额外的方程的，这里通过 correct() 函数就可以完成求解

    runTime.write();

    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;
  }
  return 0;
}
