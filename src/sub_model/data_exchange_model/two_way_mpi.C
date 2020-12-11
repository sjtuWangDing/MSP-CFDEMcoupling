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
  Foam::twoWayMPI
\*---------------------------------------------------------------------------*/

#include "./two_way_mpi.h"
#include "sub_model/liggghts_command_model/liggghts_command_model.h"

namespace Foam {

cfdemDefineTypeName(twoWayMPI)

cfdemCreateNewFunctionAdder(dataExchangeModel, twoWayMPI)

twoWayMPI::twoWayMPI(cfdemCloud& cloud)
  : dataExchangeModel(cloud),
    lmp_(nullptr),
    subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")) {

  Info << "Construct twoWayMPI..." << endl;
  MPI_Comm_dup(MPI_COMM_WORLD, &lgs_comm_);
  lmp_ = new LAMMPS_NS::LAMMPS(0, nullptr, lgs_comm_);

  // open LIGGGHTS input script
  Info << "Starting up LIGGGHTS for first time execution..."<<endl;
  Info << "Output LIGGGHTS input script:" << endl;
  const fileName lgsPath(subPropsDict_.lookup("liggghtsPath"));
  lmp_->input->file(lgsPath.c_str());
  Info << "First time LIGGGHTS execution - done" << endl;

  // get DEM time step size
  DEMts_ = lmp_->update->dt;
  Info << "First get DEM time step size: " << DEMts_ << endl;
  checkTimeStepSize();
}

twoWayMPI::~twoWayMPI() {}

//! \return 当前耦合时间步中颗粒的数量
int twoWayMPI::couple() {
  Info << "dataExchangeModel " << typeName() << ": Starting up LIGGGHTS..." << endl;
  for (const std::shared_ptr<liggghtsCommandModel>& model: cloud_.liggghtsCommandModels()) {
    liggghtsCommandModel* lgs = model.get();
    // 在 runCommand 中会判断当前的 couplingStep_ 是否满足时间步的要求
    if (lgs->runCommand(couplingStep_)) {
      int commandLines = lgs->commandLines();
      for (int j = 0; j < commandLines; ++j) {
        std::string command = lgs->getCommand(j);
        Info << "Executing command: " << command << endl;
        lmp_->input->one(command.c_str());
      }
    }
  }
  Info << "Run LIGGGHTS - done" << endl;
  // 返回颗粒数量
  return liggghts_get_maxtag(lmp_);
}

} // namespace Foam
