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

// 由于 pisoFoam 可以处理非牛顿流体，所以在此处引入了 laminarTransport
singlePhaseTransportModel laminarTransport(U, phi);

// 湍流模型
autoPtr<incompressible::turbulenceModel> turbulence
(
  incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
