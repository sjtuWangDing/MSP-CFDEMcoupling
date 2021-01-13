Info << "Reading field p\n" << endl;
volScalarField p
(
  IOobject
  (
    "p",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

Info << "Reading field U\n" << endl;
volVectorField U
(
  IOobject
  (
    "U",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

Info << "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
  IOobject
  (
    "phi",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  fvc::flux(U)
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, piso.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
  incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
