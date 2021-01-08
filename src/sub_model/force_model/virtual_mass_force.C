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

#include "./virtual_mass_force.h"

namespace Foam {

cfdemDefineTypeName(virtualMassForce);

cfdemCreateNewFunctionAdder(forceModel, virtualMassForce);

/*!
 * \brief Constructor
 * \note The initialization list should be in the same order as the variable declaration
 */
virtualMassForce::virtualMassForce(cfdemCloud& cloud)
    : forceModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      velFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("velFieldName", "U").c_str()),
      voidFractionFieldName_(
          subPropsDict_.lookupOrDefault<Foam::word>("voidFractionFieldName", "voidFraction").c_str()),
      phiFieldName_(subPropsDict_.lookupOrDefault<Foam::word>("phiFieldName", "phi").c_str()),
      U_(cloud.mesh().lookupObject<volVectorField>(velFieldName_)),
      voidFraction_(cloud.mesh().lookupObject<volScalarField>(voidFractionFieldName_)),
      phi_(cloud.mesh().lookupObject<surfaceScalarField>(phiFieldName_)) {
  createForceSubModels(subPropsDict_, kUnResolved);
}

virtualMassForce::~virtualMassForce() {}

void virtualMassForce::setForce() {
  Info << "Setting virtual mass force..." << endl;
  const volScalarField& rhoField = forceSubModel_->rhoField();
  // 计算 ddtU
  volVectorField ddtUField = fvc::ddt(U_) + fvc::div(phi_, U_);
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // 虚拟质量力
    Foam::vector virtualMassForce = Foam::vector::zero;
    if (cloud_.checkMiddleParticle(index)) {
      setMiddleParticleForceKernel(virtualMassForce, index, rhoField, ddtUField);
    }
    // 更新颗粒速度，所有处理器都必须更新，否则当颗粒从一个计算域运动到另一个计算域，新的计算域必须有颗粒的 prevU
    if (prevParticleVelMap_.end() == prevParticleVelMap_.find(index)) {
      prevParticleVelMap_.insert(std::make_pair(index, cloud_.getVelocity(index)));
    } else {
      prevParticleVelMap_[index] = cloud_.getVelocity(index);
    }
    forceSubModel_->partToArray(index, virtualMassForce, Foam::vector::zero, Foam::vector::zero, 0);
  }
  Info << "Setting virtual mass force - done" << endl;
}

void virtualMassForce::setMiddleParticleForceKernel(Foam::vector& virtualMassForce, const int index,
                                                    const volScalarField& rhoField, const volVectorField& ddtUField) {
  int findCellID = cloud_.findCellIDs()[index];
  if (findCellID >= 0) {
    std::unordered_set<int> set;                           // 颗粒扩展网格的集合
    double dt = cloud_.mesh().time().deltaT().value();     // 时间步长
    double radius = cloud_.getRadius(index);               // 颗粒半径
    double diameter = radius;                              // 颗粒直径
    double rho = rhoField[findCellID];                     // 流体密度
    double vf = 0.0;                                       // 流体空隙率
    double Ac = 0.0;                                       // Ac 系数
    Foam::vector particlePos = cloud_.getPosition(index);  // 颗粒中心坐标
    Foam::vector Up = cloud_.getVelocity(index);           // 颗粒速度
    Foam::vector prevUp = Foam::vector::zero;              // 上一个时间步颗粒速度
    Foam::vector UFluid = Foam::vector::zero;              // 背景流体速度
    Foam::vector ddtU = Foam::vector::zero;                // 背景流体加速度
    Foam::vector ddtUp = Foam::vector::zero;               // 颗粒加速度
    Foam::vector Ur = Foam::vector::zero;                  // 相对速度
    Foam::vector ddtUr = Foam::vector::zero;               // 相对加速度
    // 获取当前颗粒的 Expanded Cell 集合
    cloud_.voidFractionM().buildExpandedCellSet(set, findCellID, particlePos, radius, 6);
    // 计算背景流体速度与背景流体加速度
    double sumCore = 0.0;
    double core = 0.0;
    double cellV = 0.0;
    double sumPV = 0.0;
    double sumCV = 0.0;
    Foam::vector cellPos = Foam::vector::zero;
    for (int cellID : set) {
      if (cellID >= 0) {  // cell found
        cellPos = cloud_.mesh().C()[cellID];
        cellV = cloud_.mesh().V()[cellID];
        // 计算高斯核
        core = GaussCore(particlePos, cellPos, radius, 6);
        // 计算累计速度
        UFluid += voidFraction_[cellID] * U_[cellID] * core * cellV;
        // 计算累计加速度
        ddtU += voidFraction_[cellID] * ddtUField[cellID] * core * cellV;
        // 计算累加因数
        sumCore += voidFraction_[cellID] * core * cellV;
        // 计算累加流体体积
        sumPV += voidFraction_[cellID] * cellV;
        // 计算累加网格体积
        sumCV += cellV;
      }
    }
    UFluid /= sumCore;
    ddtU /= sumCore;
    vf = sumPV / sumCV;
    // 计算颗粒加速度 ddtUp
    if (prevParticleVelMap_.end() == prevParticleVelMap_.find(index)) {
      ddtUp = Up / dt;
    } else {
      prevUp = prevParticleVelMap_[index];
      ddtUp = (Up - prevUp) / dt;
    }
    // 计算 Ac 系数
    Ur = UFluid - Up;
    ddtUr = ddtU - ddtUp;
    Ac = sqr(mag(Ur)) / (diameter * mag(ddtUr));
    // 计算虚拟质量力
    virtualMassForce = (2.0 / 3.0) * M_PI * (2.1 - 0.132 / (0.12 + Ac * Ac)) * pow(diameter, 3) * rho * ddtUr;
    if ("B" == cloud_.modelType()) {
      virtualMassForce /= vf;
    }
    if (forceSubModel_->verbose()) {
      Pout << "vf = " << vf << endl;
      Pout << "ddtU = " << ddtU << endl;
      Pout << "ddtUp = " << ddtUp << endl;
      Pout << "ddtUr = " << ddtUr << endl;
      Pout << "Ac = " << Ac << endl;
      Pout << "virtual mass force = " << virtualMassForce << endl;
    }
  }
}

}  // namespace Foam
