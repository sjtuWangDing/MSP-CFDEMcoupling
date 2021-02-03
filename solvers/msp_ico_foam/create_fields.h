Info << "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
  IOobject
  (
    "transportProperties",
    runTime.constant(), // 从 constant 文件夹下读取 transportProperties 文件
    mesh,
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
  )
);

// 从 transportProperties 中读取有量纲的标量 nu
dimensionedScalar nu
(
  "nu",
  dimViscosity,
  transportProperties.lookup("nu")
);

Info << "Reading field p\n" << endl;
volScalarField p
(
  IOobject
  (
    "p",
    runTime.timeName(), // 从时间文件夹下读取 p
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

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
// 在 fvSolution 文件中寻找关键字 PISO，并读取 pRefCell 和 pRefValue 的值
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
