#include "fvCFD.H"
#include "fvOptions.H"
#include "pisoControl.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"

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
      fvm::ddt(U) + fvm::div(phi, U)
    + MRF.DDt(U)
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

      laminarTransport.correct();
      turbulence->correct();

      runTime.write();

      Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
          << "  ClockTime = " << runTime.elapsedClockTime() << " s"
          << nl << endl;
  }

  return 0;
}
