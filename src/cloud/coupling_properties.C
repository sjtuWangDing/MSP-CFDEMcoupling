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
  Foam::CouplingProperties
\*---------------------------------------------------------------------------*/

#include "cloud/coupling_properties.h"

const std::string CFDEM_VERSION = "cfdem-3.8.1";
const std::string LIGGGHTS_VERSION = "3.8.0";

namespace Foam {

CouplingProperties::CouplingProperties(const fvMesh& mesh, const IOdictionary& couplingPropertiesDict,
                                       const IOdictionary& liggghtsCommandsDict)
    : mesh_(mesh),
      couplingPropertiesDict_(couplingPropertiesDict),
      liggghtsCommandsDict_(liggghtsCommandsDict),
      fineParticleRatio_(0.0),
      coarseParticleRatio_(0.0),
      expandedCellScale_(0.0),
      verbose_(couplingPropertiesDict.lookupOrDefault<bool>("verbose", false)),
      solveFlow_(couplingPropertiesDict.lookupOrDefault<bool>("solveFlow", true)),
      modelType_(couplingPropertiesDict.lookupOrDefault<Foam::word>("modelType", "none").c_str()),
      turbulenceModelType_(couplingPropertiesDict.lookupOrDefault<Foam::word>("turbulenceModelType", "none").c_str()),
      allowUseSubCFDTimeStep_(couplingPropertiesDict.lookupOrDefault<bool>("allowUseSubCFDTimeStep", false)),
      couplingInterval_(couplingPropertiesDict.lookupOrDefault<int>("couplingInterval", 0)),
      checkPeriodicCells_(couplingPropertiesDict.lookupOrDefault<bool>("checkPeriodicCells", false)),
      periodicCheckRange_(Foam::vector(1, 1, 1)),
      refineMeshSkin_(couplingPropertiesDict.lookupOrDefault<double>("refineMeshSkin", 1.8)),
      refineMeshKeepInterval_(couplingPropertiesDict.lookupOrDefault<int>("refineMeshKeepInterval", 0)),
      ddtVoidFractionType_(couplingPropertiesDict.lookupOrDefault<Foam::word>("ddtVoidFractionType", "off").c_str()) {
  Info << "CFDEM coupling version: " << CFDEM_VERSION << endl;
  Info << "LIGGGHTS version: " << LIGGGHTS_VERSION << endl;

  if (couplingPropertiesDict.found("forceModels")) {
    Foam::wordList fList = couplingPropertiesDict.lookup("forceModels");
    for (Foam::wordList::iterator it = fList.begin(); it != fList.end(); ++it) {
      forceModelList_.emplace_back(it->c_str());
    }
  }

  if (couplingPropertiesDict.found("momCoupleModels")) {
    Foam::wordList mList = couplingPropertiesDict.lookup("momCoupleModels");
    for (Foam::wordList::iterator it = mList.begin(); it != mList.end(); ++it) {
      momCoupleModelList_.emplace_back(it->c_str());
    }
  }

  if (liggghtsCommandsDict.found("liggghtsCommandModels")) {
    Foam::wordList lForce = liggghtsCommandsDict.lookup("liggghtsCommandModels");
    for (Foam::wordList::iterator it = lForce.begin(); it != lForce.end(); ++it) {
      liggghtsCommandModelList_.emplace_back(it->c_str());
    }
  }

  if (refineMeshSkin_ < 1.0) {
    FatalError << "refineMeshSkin should be >= 1.0 but get " << refineMeshSkin_ << abort(FatalError);
  }

  if (refineMeshKeepInterval_ < 0) {
    FatalError << "refineMeshKeepInterval should be >= 0 but get " << refineMeshKeepInterval_ << abort(FatalError);
  }

  if ("a" != ddtVoidFractionType_ && "b" != ddtVoidFractionType_ && "off" != ddtVoidFractionType_) {
    FatalError << "Model " << ddtVoidFractionType_
               << " is not a valid choice for ddt(voidfraction). Choose a or b or off." << abort(FatalError);
  }

  fineParticleRatio_ = couplingPropertiesDict.lookupOrDefault<double>("fineParticleRatio", 3.0);

  coarseParticleRatio_ = couplingPropertiesDict.lookupOrDefault<double>("coarseParticleRatio", 0.33);

  expandedCellScale_ = couplingPropertiesDict.lookupOrDefault<double>("expandedCellScale", 6.0);

#if CFDEM_MIX_CLOUD
  fineParticleRatio_ = couplingPropertiesDict.lookupOrDefault<double>("fineParticleRatio", 3.0);
  coarseParticleRatio_ = couplingPropertiesDict.lookupOrDefault<double>("coarseParticleRatio", 0.33);
  usedForSolverIB_ = couplingPropertiesDict.lookupOrDefault<Switch>("usedForSolverIB", false);
  usedForSolverPiso_ = couplingPropertiesDict.lookupOrDefault<Switch>("usedForSolverPiso", false);
  useDynamicRefineMesh_ = couplingPropertiesDict.lookupOrDefault<Switch>("useDynamicRefineMesh", false);

  // 读取 fixed_particle_ 与来流速度
  fixedParticle_ = couplingPropertiesDict.lookupOrDefault<bool>("fixedParticle", false);
  flowVelocity_ = vector(couplingPropertiesDict.lookupOrDefault<double>("U0x", 0.0),
                         couplingPropertiesDict.lookupOrDefault<double>("U0y", 0.0),
                         couplingPropertiesDict.lookupOrDefault<double>("U0z", 0.0));

  if (fixedParticle_) {
    Info << "Using fixed particle, reading flow_velocity_: " << flowVelocity_ << endl;
  }

  if (mag(flowVelocity_) < SMALL) {
    FatalError << "mag(flowVelocity_) is zero" << abort(FatalError);
  }
#endif  // CFDEM_MIX_CLOUD

  Info << "Reading coupling properties from dictionary - done" << endl;
}

}  // namespace Foam
