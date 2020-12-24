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
  runLiggghts
\*---------------------------------------------------------------------------*/

#include "./run_liggghts.h"

namespace Foam {

cfdemDefineTypeName(runLiggghts);

cfdemCreateNewFunctionAdder(liggghtsCommandModel, runLiggghts);

//! \brief Constructor
runLiggghts::runLiggghts(cfdemCloud& cloud) : liggghtsCommandModel(cloud), subPropsDict_(cloud.liggghtsCommandsDict()) {
  std::string dictName = typeName_ + "Props";
  if (cloud.liggghtsCommandsDict().found(dictName)) {
    subPropsDict_ = cloud.liggghtsCommandsDict().subDict(dictName);
    verbose_ = subPropsDict_.lookupOrDefault<Switch>("verbose", false);
  }
  // check some switches: runFirst, runLast, runEveryCouplingStep and runEveryWriteStep
  checkTimeMode(subPropsDict_);
  checkTimeSettings(subPropsDict_);
  command_ = createCommand(baseCommand_);
}

//! \brief Destructor
runLiggghts::~runLiggghts() {}

std::string runLiggghts::getCommand(int index) const {
  return command_;
}

std::string runLiggghts::createCommand(const std::string& cmd, int interval) {
  return cmd + " " + std::to_string(interval);
}

bool runLiggghts::runCommand(int couplingStep) {
  command_ = createCommand(baseCommand_, cloud_.couplingInterval());
  return runThisCommand(couplingStep);
}

}  // namespace Foam
