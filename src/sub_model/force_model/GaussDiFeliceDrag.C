/*---------------------------------------------------------------------------*\
  CFDEMcoupling - Open Source CFD-DEM coupling

  CFDEMcoupling is part of the CFDEMproject
  www.cfdem.com
                              Christoph Goniva, christoph.goniva@cfdem.com
                              Copyright 2009-2012 JKU Linz
                              Copyright 2012-     DCS Computing GmbH, Linz
------------------------------------------------------------------------------
License
  This file is part of CFDEMcoupling.

  CFDEMcoupling is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 3 of the License, or (at your
  option) any later version.

  CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with CFDEMcoupling; if not, write to the Free Software Foundation,
  Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
\*---------------------------------------------------------------------------*/

#include "./GaussDiFeliceDrag.h"

namespace Foam {

cfdemDefineTypeName(GaussDiFeliceDrag);

cfdemCreateNewFunctionAdder(forceModel, GaussDiFeliceDrag);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
GaussDiFeliceDrag::GaussDiFeliceDrag(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      velFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("velFieldName", "U").c_str()),
      voidFractionFieldName_(
          subPropsDict_.lookupOrDefault<Foam::word>("voidFractionFieldName", "voidFraction").c_str()),
      U_(cloud.mesh().lookupObject<volVectorField>(velFieldName_)),
      voidFraction_(cloud.mesh().lookupObject<volScalarField>(voidFractionFieldName_)) {
  createForceSubModels(subPropsDict_, kUnResolved);
  CHECK(false == forceSubModel_->interpolation()) << ": GaussDiFeliceDrag model request interpolation == false";
}

GaussDiFeliceDrag::~GaussDiFeliceDrag() {}

void GaussDiFeliceDrag::setForce() {
  Info << "Setting GaussDiFeliceDrag force..." << endl;
  const volScalarField& nuField = forceSubModel_->nuField();
  const volScalarField& rhoField = forceSubModel_->rhoField();
  int findCellID = -1;           // 颗粒中心所在网格的索引
  double radius = 0.0;           // 颗粒半径
  double diameter = 0.0;         // 颗粒直径
  double nuf = 0.0;              // 流体动力粘度
  double rho = 0.0;              // 流体密度
  double dragCoefficient = 0.0;  // 颗粒阻力系数
  double vf = 0.0;               // 颗粒中心所在网格的空隙率(可以指定是否使用插值模型计算)
  double pRe = 0.0;              // 颗粒雷诺数
  double Cd = 0.0;               // 流体阻力系数 Cd = sqr(0.63 + 4.8 / sqrt(pRe))
  double Xi = 0.0;               // 模型阻力系数 Xi = 3.7 - 0.65 * exp(-sqr(1.5 - log10(pRe)) / 2)
  double magUr = 0.0;            // 相对速度值
  Foam::vector drag;             // 总阻力 = dragCoefficient * Ur
  Foam::vector Up;               // 颗粒速度
  Foam::vector Ufluid;           // 颗粒中心处流体速度(可以指定是否使用插值模型计算)
  Foam::vector Ur;               // 相对速度
  Foam::vector particlePos;      // 颗粒中心
  std::unordered_set<int> set;   // 颗粒扩展网格的集合

  // voidFractionInterpolator_.reset(
  //     interpolation<Foam::scalar>::New(
  //         subPropsDict_.lookupOrDefault("voidfractionInterpolationType", word("cellPoint")), voidFraction_)
  //         .ptr());

  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // init before calculate force.
    set.clear();
    drag = Foam::vector::zero;
    Ufluid = Foam::vector::zero;
    vf = 0.0;
    dragCoefficient = 0.0;
    findCellID = cloud_.findCellIDs()[index];
    if (findCellID >= 0) {
      particlePos = cloud_.getPosition(index);
      radius = cloud_.getRadius(index);
      // 获取当前颗粒的 Expanded Cell 集合
      cloud_.voidFractionM().buildExpandedCellSet(set, findCellID, particlePos, radius, 6);
      double sumCore = 0.0;
      double core = 0.0;
      double cellV = 0.0;
      double sumPV = 0.0;
      double sumCV = 0.0;
      Foam::vector cellPos = Foam::vector::zero;
      for (int cellID : set) {
        // for (int subCell = 0; subCell < cloud_.particleOverMeshNumber()[index]; ++subCell) {
        //   int cellID = cloud_.cellIDs()[index][subCell];
        if (cellID >= 0) {  // cell found
          cellPos = cloud_.mesh().C()[cellID];
          cellV = cloud_.mesh().V()[cellID];
          // 计算高斯核
          core = GaussCore(particlePos, cellPos, radius, 6);
          // 计算累计速度
          Ufluid += voidFraction_[cellID] * U_[cellID] * core * cellV;
          // 计算累加因数
          sumCore += voidFraction_[cellID] * core * cellV;
          // 计算累加流体体积
          sumPV += voidFraction_[cellID] * cellV;
          // 计算累加网格体积
          sumCV += cellV;
        }
      }
      // 计算平均流体速度
      Ufluid /= sumCore;
      // 计算空隙率
      vf = sumPV / sumCV;
      Up = cloud_.getVelocity(index);
      diameter = 2 * cloud_.getRadius(index);
      Ur = Ufluid - Up;
      magUr = mag(Ur);
      nuf = nuField[findCellID];
      rho = rhoField[findCellID];
      pRe = 0.0;
      Cd = 0.0;
      dragCoefficient = 0.0;
      if (magUr > 0) {
#if 0
        // 计算颗粒雷诺数
        pRe = diameter * vf * magUr / (nuf + Foam::SMALL);
        // 计算流体阻力系数
        Cd = sqr(0.63 + 4.8 / sqrt(pRe));
        // 计算模型阻力系数
        Xi = 3.7 - 0.65 * exp(-sqr(1.5 - log10(pRe)) / 2.0);
        // 计算颗粒阻力系数
        dragCoefficient = 0.125 * Cd * rho * M_PI * diameter * diameter * pow(vf, (2 - Xi)) * magUr;
#elif 1
        // 这个模型在两个 test 中，结果与实验值比较吻合
        pRe = diameter * magUr * vf / (nuf + Foam::SMALL);
        Cd = 24 * pow(9.06 / sqrt(pRe) + 1, 2) / (9.06 * 9.06);
        Xi = 3.7 - 0.65 * exp(-sqr(1.5 - log10(pRe)) / 2.0);
        dragCoefficient = 0.125 * Cd * rho * M_PI * diameter * diameter * pow(vf, (2 - Xi)) * magUr;
#elif 0
        // 计算颗粒雷诺数
        pRe = diameter * vf * magUr / (nuf + Foam::SMALL);
        if (vf > 0.8) {
          // 计算流体阻力系数
          Cd = pRe >= 1000 ? 0.44 : 24.0 * (1 + 0.15 * pow(pRe, 0.687)) / pRe;
          dragCoefficient = 0.75 * rho * vf * Cd * magUr / (diameter * Foam::pow(vf, 2.65));
        } else {
          dragCoefficient = 150 * (1 - vf) * nuf * rho / (vf * diameter * diameter) + 1.75 * magUr * rho / diameter;
        }
        dragCoefficient *= cloud_.voidFractionM().pV(radius);
#elif 0
        // 计算颗粒雷诺数
        pRe = diameter * vf * magUr / (nuf + Foam::SMALL);
        if (vf > 0.74) {
          double Wd = 0.0;
          if (vf <= 0.82) {
            Wd = 0.0214 / (4 * (vf - 0.7463) + 0.0044) - 0.576;
          } else if (vf > 0.97) {
            Wd = 32.8295 * vf - 31.8295;
          } else {
            Wd = 0.0038 / (4 * sqr(vf - 0.7789) + 0.004) - 0.0101;
          }
          Cd = pRe >= 1000 ? 0.44 : 24.0 * (1 + 0.15 * pow(pRe, 0.687)) / pRe;
          dragCoefficient = 0.75 * rho * vf * magUr * Cd * Wd / diameter;
        } else {
          dragCoefficient = 150 * (1 - vf) * nuf * rho / (vf * diameter * diameter) + 1.75 * magUr * rho / diameter;
        }
        dragCoefficient *= cloud_.voidFractionM().pV(radius);
#endif
        if ("B" == cloud_.modelType()) {
          dragCoefficient /= vf;
        }
        // 计算总阻力
        drag = dragCoefficient * Ur;
      }
      if (forceSubModel_->verbose()) {
        Pout << "index = " << index << endl;
        Pout << "findCellID = " << findCellID << endl;
        Pout << "Up = " << Up << endl;
        Pout << "Ur = " << Ur << endl;
        Pout << "diameter = " << diameter << endl;
        Pout << "rho = " << rho << endl;
        Pout << "nuf = " << nuf << endl;
        Pout << "voidFraction = " << vf << endl;
        Pout << "pRe = " << pRe << endl;
        Pout << "Xi = " << Xi << ", " << pow(vf, (2 - Xi)) << endl;
        Pout << "Cd = " << Cd << endl;
        Pout << "dragCoefficient = " << dragCoefficient << endl;
        Pout << "drag (total) = " << drag << endl;
      }
    }
    base::MPI_Barrier();
    // write particle data to global array
    forceSubModel_->partToArray(index, drag, Foam::vector::zero, Ufluid, dragCoefficient);
  }
  Info << "Setting GaussDiFeliceDrag - done" << endl;
}

}  // namespace Foam
