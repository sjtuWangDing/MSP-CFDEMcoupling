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

Info << "Reading density field rho\n" << endl;
volScalarField rho(
  IOobject(
    "rho",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
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

Info << "Reading voidFraction field voidFraction (= Vgas / Vparticle) of fine and middle particles\n" << endl;
volScalarField voidFraction(
  IOobject(
    "voidFraction",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1.0)  // 如果在时间文件夹下没有 voidFraction，那么其应该被初始化为 1.0
);

Info << "Reading volumeFraction field volumeFraction (= Vgas / Vparticle) of coarse particles\n" << endl;
volScalarField volumeFraction(
  IOobject(
    "volumeFraction",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1.0)  // 如果在时间文件夹下没有 volumeFraction，那么其应该被初始化为 1.0
);

Info << "Reading momentum exchange field Ksl\n" << endl;
volScalarField Ksl(
  IOobject(
    "Ksl",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.0)  // F / (V * Ur) ~ kg * m^-3 * s^-1
);

Info << "Reading particle velocity field Us\n" << endl;
volVectorField Us(
  IOobject(
    "Us",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0), Foam::vector::zero)
);

Info << "Reading field interface for dynamic refine mesh\n" << endl;
volScalarField interface(
  IOobject(
    "interface",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

Info << "Reading field phiIB\n" << endl;
volScalarField phiIB(
  IOobject(
    "phiIB",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
);

Info << "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi(
  IOobject(
    "phi",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  linearInterpolate(U * voidFraction) & mesh.Sf()
);

Info << "Reading/calculating face flux field phiByVoidFraction\n" << endl;
surfaceScalarField phiByVoidFraction(
  IOobject(
    "phiByVoidFraction",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  linearInterpolate(U) & mesh.Sf()
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence(
  incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
