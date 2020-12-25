Info << "Reading field p\n" << endl;
volScalarField p(
  IOobject(
    "p",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

Info << "Reading physical velocity field U\n" << endl;
volVectorField U(
  IOobject(
    "U",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

Info << "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi(
  IOobject(
    "phi",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  fvc::flux(U)
);

Info << "Reading density field rho\n" << endl;
volScalarField rho(
  IOobject(
    "rho",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh
);

Info << "Reading gravity field g\n" << endl;
uniformDimensionedVectorField g(
  IOobject
  (
    "g",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
  )
);

Info << "Reading volumeFraction field volumeFraction = (Vgas / Vparticle)\n" << endl;
volScalarField volumeFraction(
  IOobject(
    "volumeFraction",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1.0)
);

Info << "Reading field interface\n" << endl;
volScalarField interface(
  IOobject(
    "interface",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

Info << "Reading field refineMeshKeepStep\n" << endl;
volScalarField refineMeshKeepStep(
  IOobject(
    "refineMeshKeepStep",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

Info << "Reading field phiIB\n" << endl;
volScalarField phiIB(
  IOobject(
    "phiIB",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence(
  incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
