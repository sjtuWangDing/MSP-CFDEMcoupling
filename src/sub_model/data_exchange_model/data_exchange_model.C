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

Description
  This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
  and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

Class
  dataExchangeModel
\*---------------------------------------------------------------------------*/

#include "./data_exchange_model.h"

namespace Foam {

cfdemDefineTypeName(dataExchangeModel)

cfdemDefineNewFunctionMap(dataExchangeModel)

cfdemDefineConstructNewFunctionMap(dataExchangeModel)

cfdemDefineDestroyNewFunctionMap(dataExchangeModel)

cfdmeDefineBaseTypeNew(autoPtr, dataExchangeModel, (cfdemCloud& cloud, const dictionary& dict), dict, (cloud))

//! \brief Constructor
dataExchangeModel::dataExchangeModel(cfdemCloud& cloud)
  : cloud_(cloud),
    // 在初始化 dataExchangeModel 的时候，记录下当前的流体时间步为 time index
    timeIndexOffset_(cloud.mesh().time().timeIndex()),
    // 初始化耦合时间步
    couplingStep_(0),
    // 初始化 DEM 时间步长，在具体的模型中读入具体的值
    // twoWayMPI: 通过 LAMMP 读入
    // twoWayFile: 通过字典文件读入
    DEMts_(-1.0) {
  Info << "dataExchangeMode: timeIndexOffset_ = " << timeIndexOffset_ << endl;
  Info << "dataExchangeMode: couplingStep_ = " << couplingStep_ << endl;
}

//! \brief Destructor
dataExchangeModel::~dataExchangeModel() {}

/*!
 * \brief 检查时间步长是否满足耦合要求（在初始化模型的时候调用）
 * \note (1) 耦合时间步长应该大于 or 等于 CFD 时间步长
 *       (2) 耦合时间步长应该为 CFD 时间步长的整数倍
 *       (3) 如果耦合时间步长 != CFD 时间步长，则要求 allowUseSubCFDTimeStep() 为 true
 */
void dataExchangeModel::checkTimeStepSize() const {
  scalar CFDts = cloud_.mesh().time().deltaT().value(); // 获取 CFD 时间步长
  Info << "CFD time step size: " << CFDts << endl;
  // 耦合时间步长 >= CFD 时间步长
  if (CFDts > couplingTime() + Foam::SMALL) {
    Info << "CFD time step size = " << CFDts << endl;
    Info << "DEM time step size = " << DEMts_ << endl;
    Info << "Coupling interval = " << cloud_.cProps().couplingInterval() << endl;
    FatalError << "CFD time-step bigger than coupling time (= DEM time step * coupling interval)!\n"
      << abort(FatalError);
  }
  // 耦合时间步长应该为 CFD 时间步长的整数倍
  if (std::fabs((std::round(couplingTime() / CFDts) * CFDts) - couplingTime()) > Foam::SMALL) {
    Info << "CFD time step size = " << CFDts << endl;
    Info << "DEM time step size = " << DEMts_ << endl;
    Info << "Coupling interval = " << cloud_.cProps().couplingInterval() << endl;
    FatalError << "Coupling time (= DEM time step * coupling interval) is not a multiple of  CFD time-step!\n"<< abort(FatalError);
  }
  // 如果耦合时间步长 != CFD 时间步长
  if (couplingTime() > CFDts + Foam::SMALL) {
    // 如果不允许使用 CFD sub time step 则报错
    if (!cloud_.cProps().allowUseSubCFDTimeStep()) {
      FatalError << "Your models require: CFD time-step = coupling interval (= DEM time step * coupling interval)! \n" << abort(FatalError);
    }
    Warning << "You are using sub-time-steps (i.e. CFD time-step < coupling time)" << endl;
  }
}

/*!
 * \brief 因为耦合时间步长 = 流体时间步长的整数倍，所以 timeStepFraction() 用于计算每个流体时间步在耦合时间步中的所占比例，
 *        如果 couplingTime() == 3 * CFDts，那么每一个耦合时间步由 3 个流体时间步构成，
 *        那么这三个流体时间步的 timeStepFraction() 分别返回 0, 0.333333, 0.666666
 */
double dataExchangeModel::timeStepFraction() const {
  Info << "timeIndex: " << cloud_.mesh().time().timeIndex() << ", couplingStep: " << couplingStep_ << endl;
  // 如果当前流体时间步是耦合时间步中第一个流体时间步
  if (checkValidCouplingStep()) {
    return 0;
  }
  scalar CFDts = cloud_.mesh().time().deltaT().value();
  label timeIndex = cloud_.mesh().time().timeIndex();
  // (timeIndex - timeIndexOffset_ - 1) * CFDts: 当前流体时间步起始时间，Eg: 0.00, 0.01, 0.02, 0.03
  // (couplingStep_ - 1) * couplingTime(): 耦合时间步截止到 (couplingStep_ - 1) 的计算时间
  return ((timeIndex - timeIndexOffset_ - 1) * CFDts - (couplingStep_ - 1) * couplingTime()) / couplingTime();
}

/*!
 * \brief 因为耦合时间步长 = 流体时间步长的整数倍，而耦合发生在耦合时间步中的第一个流体时间步中，
 *        所以判断当前流体时间步是否同时也是耦合时间步
 */
bool dataExchangeModel::checkValidCouplingStep() const {
  scalar CFDts = cloud_.mesh().time().deltaT().value();
  // OpenFOAM 中 timeIndex 从 1 开始
  label timeIndex = cloud_.mesh().time().timeIndex();
  // (timeIndex - timeIndexOffset_ - 1) * CFDts: 当前流体时间步起始时间，Eg: 0.00, 0.01, 0.02, 0.03
  // couplingStep_ * couplingTime(): 当前耦合时间步的起始时间，Eg: 0.00, 0.03, 0.06
  if ((timeIndex - timeIndexOffset_ - 1) * CFDts + Foam::SMALL > couplingStep_ * couplingTime()) {
    return true;
  }
  return false;
}

} // namespace Foam
