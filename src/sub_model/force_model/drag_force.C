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

#include "./drag_force.h"

namespace Foam {

cfdemDefineTypeName(dragForce);

cfdemCreateNewFunctionAdder(forceModel, dragForce);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
dragForce::dragForce(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      dragModelName_(subPropsDict_.lookupOrDefault<Foam::word>("dragModelName", "DiFelice").c_str()),
      velFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("velFieldName", "U").c_str()),
      UsFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("granVelFieldName", "Us").c_str()),
      voidFractionFieldName_(
          subPropsDict_.lookupOrDefault<Foam::word>("voidFractionFieldName", "voidFraction").c_str()),
      U_(cloud.mesh().lookupObject<volVectorField>(velFieldName_)),
      UsField_(cloud.mesh().lookupObject<volVectorField>(UsFieldName_)),
      voidFraction_(cloud.mesh().lookupObject<volScalarField>(voidFractionFieldName_)) {
  createForceSubModels(subPropsDict_, kUnResolved);
  size_t dragModelHashValue = strHasher_(dragModelName_);
  if (DiFeliceHashValue_ != dragModelHashValue && AbrahamHashValue_ != dragModelHashValue &&
      SchillerNaumannHashValue_ != dragModelHashValue && GidaspowHashValue_ != dragModelHashValue &&
      SyamlalObrienHashValue_ != dragModelHashValue) {
    FatalError << __func__ << ": wrong drag model name: " << dragModelName_ << abort(FatalError);
  }
  Info << __func__ << ": choose " << dragModelName_ << " drag force model." << endl;
}  // namespace Foam

dragForce::~dragForce() {}

void dragForce::setForce() {
  Info << "Setting " << dragModelName_ << " drag force..." << endl;
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

  // #include "resetVoidfractionInterpolator.H"
  // #include "resetUInterpolator.H"
  UInterpolator_.reset(
      interpolation<Foam::vector>::New(subPropsDict_.lookupOrDefault("UInterpolationType", word("cellPointFace")), U_)
          .ptr());
  voidFractionInterpolator_.reset(
      interpolation<Foam::scalar>::New(
          subPropsDict_.lookupOrDefault("voidfractionInterpolationType", word("cellPoint")), voidFraction_)
          .ptr());

  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // init
    drag = Foam::vector::zero;
    Ufluid = Foam::vector::zero;
    dragCoefficient = 0.0;
    // 这里使用 findCellID，不要使用 cellIDs，因为如果当前求解器不包含颗粒，则 cellIDs[index].mSize() == 0
    findCellID = cloud_.findCellIDs()[index];
    if (findCellID >= 0) {
      if (forceSubModel_->interpolation()) {
        // 获取颗粒中心的坐标, 将颗粒中心所在网格的空隙率和流体速度插值到颗粒中心处
        Foam::vector pos = cloud_.getPosition(index);
        vf = voidFractionInterpolator_().interpolate(pos, findCellID);
        Ufluid = UInterpolator_().interpolate(pos, findCellID);
        // 确保插值后颗粒中心的空隙率有意义
        vf = vf > 1.0 ? 1.0 : vf;
        vf = vf < Foam::SMALL ? Foam::SMALL : vf;
      } else {
        // 不使用插值模型，直接设置颗粒中心处的值为颗粒中心所在网格的值
        vf = voidFraction_[findCellID];
        Ufluid = U_[findCellID];
      }
      // 初始化
      Up = cloud_.getVelocity(index);
      radius = cloud_.getRadius(index);
      diameter = 2 * radius;
      Ur = Ufluid - Up;
      magUr = mag(Ur);
      nuf = nuField[findCellID];
      rho = rhoField[findCellID];
      pRe = 0.0;
      Cd = 0.0;
      dragCoefficient = 0.0;
      if (magUr > 0) {
        size_t dragModelHashValue = strHasher_(dragModelName_);
        // 计算 drag force
        if (DiFeliceHashValue_ == dragModelHashValue) {
          // DiFelice drag model
          pRe = diameter * vf * magUr / (nuf + Foam::SMALL);
          Cd = sqr(0.63 + 4.8 / sqrt(pRe));
          Xi = 3.7 - 0.65 * exp(-sqr(1.5 - log10(pRe)) / 2.0);
          dragCoefficient = 0.125 * Cd * rho * M_PI * diameter * diameter * pow(vf, (2 - Xi)) * magUr;
        } else if (AbrahamHashValue_ == dragModelHashValue) {
          // Abraham drag model
          pRe = diameter * magUr / (nuf + Foam::SMALL);
          Cd = 24 * pow(9.06 / sqrt(pRe) + 1, 2) / (9.06 * 9.06);
          Xi = 3.7 - 0.65 * exp(-sqr(1.5 - log10(pRe)) / 2.0);
          dragCoefficient = 0.125 * Cd * rho * M_PI * diameter * diameter * pow(vf, (2 - Xi)) * magUr;
        } else if (SchillerNaumannHashValue_ == dragModelHashValue) {
          // Schiller-Naumann drag model
          pRe = diameter * vf * magUr / (nuf + Foam::SMALL);
          Cd = pRe >= 1000 ? 0.44 : 24.0 * (1 + 0.15 * pow(pRe, 0.687)) / pRe;
          dragCoefficient = 0.125 * Cd * rho * M_PI * diameter * diameter * magUr;
        } else if (GidaspowHashValue_ == dragModelHashValue) {
          // Gidaspow drag model
          pRe = diameter * vf * magUr / (nuf + Foam::SMALL);
          if (vf > 0.8) {
            Cd = pRe >= 1000 ? 0.44 : 24.0 * (1 + 0.15 * pow(pRe, 0.687)) / pRe;
            dragCoefficient = 0.75 * rho * vf * Cd * magUr / (diameter * Foam::pow(vf, 2.65));
          } else {
            dragCoefficient = 150 * (1 - vf) * nuf * rho / (vf * diameter * diameter) + 1.75 * magUr * rho / diameter;
          }
          dragCoefficient *= cloud_.voidFractionM().pV(radius);
        } else if (SyamlalObrienHashValue_ == dragModelHashValue) {
          // Syamlal-Obrien drag model
          pRe = diameter * vf * magUr / (nuf + Foam::SMALL);
          double Vrs = 0.0;
          double A = 0.0;
          double B = 0.0;
          if (vf <= 0.85) {
            A = pow(vf, 4.14);
            B = 0.8 * pow(vf, 4.14);
          } else {
            A = pow(vf, 4.14);
            B = pow(vf, 2.65);
          }
          Vrs = 0.5 * (A - 0.06 * pRe + sqrt(sqr(0.06 * pRe) + 0.12 * pRe * (2 * B - A) + A * A));
          Cd = sqr(0.63 + 4.8 / sqrt(pRe / Vrs));
          dragCoefficient = 0.75 * vf * rho * Cd * magUr / (sqr(Vrs) * diameter);
          dragCoefficient *= cloud_.voidFractionM().pV(radius);
        }
        if ("B" == cloud_.modelType()) {
          dragCoefficient /= vf;
        }
        // 计算总阻力
        drag = dragCoefficient * Ur;
      }
      if (forceSubModel_->verbose() && 0 == index) {
        Pout << "index = " << index << endl;
        Pout << "findCellID = " << findCellID << endl;
        Pout << "Up = " << Up << endl;
        Pout << "Ur = " << Ur << endl;
        Pout << "diameter = " << diameter << endl;
        Pout << "rho = " << rho << endl;
        Pout << "nuf = " << nuf << endl;
        Pout << "voidFraction = " << vf << endl;
        Pout << "pRe = " << pRe << endl;
        Pout << "Cd = " << Cd << endl;
        Pout << "dragCoefficient = " << dragCoefficient << endl;
        Pout << "drag (total) = " << drag << endl;
      }
    }
    // write particle data to global array
    forceSubModel_->partToArray(index, drag, Foam::vector::zero, Ufluid, dragCoefficient);
  }
  Info << "Setting " << dragModelName_ << " drag force - done" << endl;
}

}  // namespace Foam
