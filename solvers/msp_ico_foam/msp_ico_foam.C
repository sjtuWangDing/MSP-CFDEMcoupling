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
  mspIcoFoam

Description
  Transient solver for incompressible, laminar flow of Newtonian fluids.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"

  pisoControl piso(mesh);

  #include "create_fields.h"
  #include "initContinuityErrs.H"

  Info << "\nStarting time loop\n" << endl;

  while (runTime.loop())
  {
    Info << "Time = " << runTime.timeName() << nl << endl;

    #include "CourantNo.H"

    // Momentum predictor
    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      - fvm::laplacian(nu, U)
    );

    if (piso.momentumPredictor())
    {
      solve(UEqn == -fvc::grad(p));
    }

    // 离散后的动量方程是一个线性方程组，其包含系数矩阵(A)和右边源项(As)
    // Aa * U = As - grad(p)
    // 将系数矩阵分解为对角线和非对角线
    // Aa = Ad + An
    // Ad * U + An * U = As - grad(p)
    // 定义 Ah
    // Ad * U = As - An * U - grad(p)
    // Ad * U = Ah - grad(p)
    // 定义 HbyA
    // HbyA = Ah / Ad
    // U = HbyA - grad(p) / Ad
    // 压力泊松方程
    // 0 = grad(U) = grad(HbyA - grad(p) / Ad)
    // grad(HbyA) = grad(grad(p) / Ad)
    // 求解压力泊松方程，得到 p
    // 根据上面 U = HbyA - grad(p) / Ad，以及新的 p 值，修正通量场 phi
    // phi = (HbyA - grad(p)) & mesh.Sf()

    // PISO loop
    while (piso.correct())
    {
      // rAU == 1 / Ad，是标量场，且已知，对应于线性方程组主对角线系数
      volScalarField rAU(1.0 / UEqn.A());

      // HbyA = Ah / Ad，是矢量场，且已知，因为 Ah = As - An * U，这里 U 是预测速度(未修正)
      // 等价于：
      // volVectorField HbyA("HbyA", U);
      // HbyA = rAU * UEqn.H();
      volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));

      // 定义 HbyA 的通量，是面标量场，且已知
      // HbyA = fvc::interpolate(HbyA) & mesh.Sf() + correction(修正项)
      surfaceScalarField phiHbyA
      (
          "phiHbyA",
          fvc::flux(HbyA) // fvc::interpolate(HbyA) & mesh.Sf()
        + fvc::interpolate(rAU) * fvc::ddtCorr(U, phi) // correction
      );

      adjustPhi(phiHbyA, U, p);

      // Update the pressure BCs to ensure flux consistency
      constrainPressure(p, U, phiHbyA, rAU);

      // Non-orthogonal pressure corrector loop
      while (piso.correctNonOrthogonal())
      {
        // p 方程
        fvScalarMatrix pEqn
        (
          // rAU 是体标量场
          fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
        );
        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

        if (piso.finalNonOrthogonalIter())
        {
          // 注意，这里修正的是通量场 phi，而不是速度，因为我们关心的是通量
          // phi = (HbyA - grad(p)) & mesh.Sf()
          // pEqn.flux() 返回的就是 grad(p) & mesh.Sf()
          phi = phiHbyA - pEqn.flux();
        }
      }
      #include "continuityErrs.H"

      // 通过新的压力场 p 更新速度场 U
      U = HbyA - rAU * fvc::grad(p);
      U.correctBoundaryConditions();

      // 我们在 piso 修正步中，每一次速度场都别更新，但是并没有更新 UEqn.H()，我们知道虽然 UEqn.H() 是速度场的函数，但是在小库朗数或小时间步的情况下，我们认为其变化不大
    }

    runTime.write();

    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;
  }
  Info << "End\n" << endl;
  return 0;
}


// ************************************************************************* //
