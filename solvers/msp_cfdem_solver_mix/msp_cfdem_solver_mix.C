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

#include "dynamicFvMesh.H"
#include "fvCFD.H"
#include "fvOptions.H"
#include "pisoControl.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"

#include "cfdem_tools/cfdem_tools.h"
#include "cloud/cfdem_cloud_mix.h"

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
  Foam::cfdemCloudMix particleCloud(mesh);
  Foam::cfdemTools::checkModelType(particleCloud);
  std::string modelType = particleCloud.modelType();

  Info << "\nStarting time loop\n" << endl;

  while (runTime.loop()) {
    Info << "\nTime = " << runTime.timeName() << nl << endl;

    #include "CourantNo.H"

    particleCloud.evolve(U, voidFraction, volumeFraction, Us, Ksl, interface);

    // 这里不需要计算 U 的通量，因为 phiByVoidFraction 在上个时间步的压力迭代中被修正，满足连续性方程
    // phiByVoidFraction = fvc::interpolate(U) & mesh.Sf();
    surfaceScalarField voidFractionFace = fvc::interpolate(voidFraction);
    // 需要重新计算 phi，因为 voidFraction 被更新了
    phi = voidFractionFace * phiByVoidFraction;

    Info << "Solver level total Eulerian momentum exchange:" << endl;
    if (modelType == "none") {
      volScalarField particleVoidFraction(voidFraction - 1.0);
      dimensionedScalar totalParticleVoidFraction = gSum(particleVoidFraction);
      Info << "  total particle voidFraction:  " << totalParticleVoidFraction.value() << endl;
    }
    dimensionedScalar totalKsl = gSum(Ksl);
    Info << "  total Ksl:  " << totalKsl.value() << endl;
    volVectorField fImp(Ksl * (Us - U));
    particleCloud.scaleWithVcell(fImp);
    dimensionedVector totalImplForce = gSum(fImp);
    Info << "  total impl force:  " << totalImplForce.value() << endl;

    if (particleCloud.solveFlow()) {
      // U prev equation
      fvVectorMatrix UEqnPre
      (
          // voidFraction * fvm::ddt(U)
          fvm::ddt(voidFraction, U) - fvm::Sp(fvc::ddt(voidFraction), U)
        + fvm::div(phi, U) - fvm::Sp(fvc::div(phi), U)
        + particleCloud.divVoidFractionTau(U, voidFraction)
        ==
        - fvm::Sp(Ksl / rho, U)
        + fvOptions(U)
      );
      UEqnPre.relax();
      fvOptions.constrain(UEqnPre);

      if (piso.momentumPredictor()) {
        if (modelType == "B" || modelType == "Bfull") {
          // modelType 为 "B" or "Bfull" 时, 压力项中不需要乘以空隙率
          solve(UEqnPre == -fvc::grad(p) + Ksl / rho * Us);
        } else if ("A" == modelType) {
          // modelType 为 "A" 时, 压力项中需要乘以空隙率
          solve(UEqnPre == -voidFraction * fvc::grad(p) + Ksl / rho * Us);
        } else if ("none" == modelType) {
          // modelType 为 "none" 时，直接求解压力项
          solve(UEqnPre == -fvc::grad(p));
        } else {
          FatalError << __func__ << ": Not implement for modelType = " << modelType << abort(FatalError);
        }
        fvOptions.correct(U);
      }

      phiByVoidFraction = fvc::flux(U);
      phi = voidFractionFace * phiByVoidFraction;
      // U equation
      fvVectorMatrix UEqn
      (
          // voidFraction * fvm::ddt(U)
          fvm::ddt(voidFraction, U) - fvm::Sp(fvc::ddt(voidFraction), U)
        + fvm::div(phi, U) - fvm::Sp(fvc::div(phi), U)
        + particleCloud.divVoidFractionTau(U, voidFraction)
        ==
        - fvm::Sp(Ksl / rho, U)
        + fvOptions(U)
      );
      UEqn.relax();
      fvOptions.constrain(UEqn);
      particleCloud.calcFictitiousForce(U, rho, volumeFraction, fictitiousForce);
      if (modelType == "B" || modelType == "Bfull") {
        // modelType 为 "B" or "Bfull" 时, 压力项中不需要乘以空隙率
        solve(UEqn == -fvc::grad(p) + Ksl / rho * Us + fictitiousForce);
      } else if ("A" == modelType) {
        // modelType 为 "A" 时, 压力项中需要乘以空隙率
        solve(UEqn == -voidFraction * fvc::grad(p) + Ksl / rho * Us + fictitiousForce);
      } else if ("none" == modelType) {
        // modelType 为 "none" 时，直接求解压力项
        solve(UEqn == -fvc::grad(p) + fictitiousForce);
      } else {
        FatalError << __func__ << ": Not implement for modelType = " << modelType << abort(FatalError);
      }
      fvOptions.correct(U);
      /*
        (1) model A
          - 离散后的动量方程是一个线性方程组，其包含系数矩阵(Aa)和右边源项(As)
            Aa * U = As - vf * grad(p) + Ksl / rho * Us

          - 将系数矩阵分解为对角线和非对角线
            Aa = Ad + An
            Ad * U + An * U = As - vf * grad(p) + Ksl / rho * Us
            Ad * U = (As - An * U) - vf * grad(p) + Ksl / rho * Us
            其中 UEqn.H() = As - An * U

          - 定义 rAU 和 HbyA
            rAU = 1 / Ad
            HbyA = rAU * UEqn.H()

          - U 表达式
            U = HbyA - rAU * vf * grad(p) + rAU * (Ksl / rho) * Us

          - 定义 HbyA 的通量 phiHbyA
            phiHbyA = HbyA & mesh().Sf()

          - 定义 Ksl / rho * Us 的通量 phiUs
            phiUs = rAU * (Ksl / rho) * Us & mesh().Sf()

          - 计算总通量 phiHbyA
            phiHbyA += phiUs

          - 推导压力泊松方程
            将 U 表达式带入连续方程: ddt(vf) + div(vf * U) = 0;
              ddt(vf) + div( vf * (HbyA + rAU * (Ksl / rho) * Us) ) == div( vf * rAU * vf * grad(p) )
            如果使用通量 phiHbyA 表示，即
              ddt(vf) + div( fvc::interpolate(vf) * phiHyA ) == div( vf * rAU * vf * grad(p) )

          - 求解压力泊松方程，得到 p

          - 将 p 带入 U = HbyA - vf * grad(p) + Ksl / rho * Us，计算新的 U，在代码中更新的是通量 phi
       */

      // PISO loop
      while (piso.correct()) {
        // - 定义 rAU，此处为定义体标量场 1 / Ad，其在压力修正步中保持不变
        volScalarField rAU(1.0 / UEqn.A());

        // - 定义 HbyA，Ad 在压力修正步中保持不变，但是 HbyA 是会更新的
        // - UEqn.H() = As - An * U
        volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));

        // - 计算 HbyA 的通量场(不包含 Ksl / rho * Us 的通量)
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            fvc::flux(HbyA)  // fvc::interpolate(HbyA) & mesh.Sf()
#if 1
          + fvc::interpolate(rAU * voidFraction) * fvc::ddtCorr(U, phiByVoidFraction)  // correction
#else
          + fvc::interpolate(rAU) * fvc::ddtCorr(U, phi)  // correction
#endif
        );

        // - 定义 Ksl / rho * Us 的通量
        surfaceScalarField phiUs("phiUs", fvc::interpolate(Us) & mesh.Sf());

        // - 将 Us 通量与 HbyA 通量相加，计算压力修正方程中的总通量
        phiHbyA += fvc::interpolate(rAU) * (fvc::interpolate(Ksl / rho) * phiUs);

        // - 定义 fictitiousForce 的通量
        surfaceScalarField phiFictitiousForce(fvc::interpolate(fictitiousForce) & mesh.Sf());

        // - 将 fictitiousForce 通量与 HbyA 通量相加，计算压力修正方程中的总通量
        phiHbyA += fvc::interpolate(rAU) * phiFictitiousForce;

        // - 在求解压力泊松方程的时候, 如果压力全部是 Neumann 边界条件(即第二类边界条件)，需要满足相容性条件，这里修正的是通量 phiHbyA
        // 在 adjustPhi 函数中，第二个参数必须使用 U，而不能使用 HbyA，因为在 adjustPhi 函数中，需要通过 U 获取 boundaryField
        // Ref: https://cfd-china.com/topic/501/%E5%85%B3%E4%BA%8Ecorrectphi-h%E8%BF%99%E4%B8%AA%E5%87%BD%E6%95%B0/9
        adjustPhi(phiHbyA, U, p);

        // Update the pressure BCs to ensure flux consistency
        if (modelType == "A") {
          volScalarField tempRAUVoidFraction = rAU * voidFraction;
          constrainPressure(p, HbyA, phiHbyA, tempRAUVoidFraction);
        } else {
          constrainPressure(p, HbyA, phiHbyA, rAU);
        }

        volScalarField rAUVoidFraction("voidFraction/AU", rAU * voidFraction);
        if (modelType == "A") {
          rAUVoidFraction = volScalarField("voidFraction/AU", rAU * voidFraction * voidFraction);
        }

        // Non-orthogonal pressure corrector loop
        // 通过迭代求解压力泊松方程, 并对方程中的拉普拉斯项进行非正交修正
        // 注意: 由于压力泊松方程中存在压力的拉普拉斯项, 当网格非正交的时候会出现显式源项,
        // 显式源项会在正交修正迭代后进行更新
        while (piso.correctNonOrthogonal()) {
          // p equation
          fvScalarMatrix pEqn
          (
              fvm::laplacian(rAUVoidFraction, p)
            ==
              fvc::div(voidFractionFace * phiHbyA)
            + particleCloud.ddtVoidFraction()
          );
          pEqn.setReference(pRefCell, pRefValue);

          pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

          if (piso.finalNonOrthogonalIter()) {
            // 更新通量场 phiByVoidFraction，这里 phiByVoidFraction 就是 U 的通量场，
            // 同时 pEqn.flux() 计算的是 p 方程中 rAUVoidFraction * grad(p) & Sf
            // 这里需要除以 voidFractionFace，因为 pEqn.flux() 修正的其实是 voidFractionFace * phiHbyA
            // Ref: http://www.dyfluid.cn/theory.pdf
            phiByVoidFraction = phiHbyA - pEqn.flux() / voidFractionFace;
          }
        }  // non-orthogonal corrector loop
        phi = voidFractionFace * phiByVoidFraction;

        #include "continuityErrorPhiPU.H"

        if (modelType == "B" || modelType == "Bfull") {
          U = HbyA - rAU * fvc::grad(p) + Ksl / rho * Us * rAU + fictitiousForce * rAU;
        } else if ("A" == modelType) {
          U = HbyA - voidFraction * rAU * fvc::grad(p) + Ksl / rho * Us * rAU + fictitiousForce * rAU;
        } else if ("none" == modelType) {
          U = HbyA - rAU * fvc::grad(p) + fictitiousForce * rAU;
        }
        U.correctBoundaryConditions();
        fvOptions.correct(U);
      }  // piso loop
      laminarTransport.correct();
      turbulence->correct();

#if 1
      Info << "mspCfdemSolverMix: after piso loop and recalculate fictitiousForce..." << endl;
      particleCloud.calcFictitiousForce(U, rho, volumeFraction, fictitiousForce);
#else
      // 通过颗粒速度以及phiIB修正速度与压力
      particleCloud.calcVelocityCorrection(p, U, phiIB, phi, voidFraction, volumeFraction);
      fvOptions.correct(U);
      dimensionedVector totalUCorrectedByPhiIB = gSum(volVectorField(fvc::grad(phiIB) / voidFraction));
      Info << "After calcVelocityCorrection, total U corrected by phiIB = " << totalUCorrectedByPhiIB.value() << nl << endl;
#endif

      // recalculate phiByVoidFraction for next step
      phiByVoidFraction = fvc::flux(U);
    }  // solve flow
    runTime.write();
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << ", ClockTime = " << runTime.elapsedClockTime() << " s" << endl;
  }  // runtime loop

  Info << "mspCfdemSolverMix - done\n" << endl;

  return 0;
}
