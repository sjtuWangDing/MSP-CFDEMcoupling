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
    pimpleFoam

Description
    Large time-step transient solver for incompressible, turbulent flow, using
    the PIMPLE (merged PISO-SIMPLE) algorithm.

    pimpleFoam 与 pisoFoam 相比较，核心在于适用于大时间步长下的求解，可以不满足
    库朗数 < 1 的条件，因为 pimpleFoam 引入了一个外循环 pimple.loop()

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // pimple foam 引入了外层循环，即在每个循环中引入了多个亚循环
        // 在每个时间步中，会求解若干次动量方程以及压力修正的过程
        while (pimple.loop())
        {
          MRF.correctBoundaryVelocity(U);

          tmp<fvVectorMatrix> tUEqn
          (
              fvm::ddt(U) + fvm::div(phi, U)
            + MRF.DDt(U)
            + turbulence->divDevReff(U)
            ==
              fvOptions(U)
          );
          fvVectorMatrix& UEqn = tUEqn.ref();

          UEqn.relax();

          fvOptions.constrain(UEqn);

          if (pimple.momentumPredictor())
          {
            solve(UEqn == -fvc::grad(p));
            fvOptions.correct(U);
          }

          // --- Pressure corrector loop
          while (pimple.correct())
          {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            MRF.makeRelative(phiHbyA);

            adjustPhi(phiHbyA, U, p);

            tmp<volScalarField> rAtU(rAU);

            if (pimple.consistent())
            {
              rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
              phiHbyA +=
                  fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
              HbyA -= (rAU - rAtU())*fvc::grad(p);
            }

            if (pimple.nCorrPISO() <= 1)
            {
              tUEqn.clear();
            }

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAtU(), MRF);

            // Non-orthogonal pressure corrector loop
            while (pimple.correctNonOrthogonal())
            {
              // Pressure corrector
              fvScalarMatrix pEqn
              (
                  fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
              );

              pEqn.setReference(pRefCell, pRefValue);

              pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

              if (pimple.finalNonOrthogonalIter())
              {
                  phi = phiHbyA - pEqn.flux();
              }
            }

            #include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            U = HbyA - rAtU()*fvc::grad(p);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
          }
            // 是否在每一个亚循环中修正非牛顿流体和湍流模型的参数
            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

      runTime.write();

      Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
          << "  ClockTime = " << runTime.elapsedClockTime() << " s"
          << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
