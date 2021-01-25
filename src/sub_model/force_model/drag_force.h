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
  The force model performs the calculation of forces (e.g. fluid-particle
  interaction forces) acting on each DEM particle.

Syntax
  forceModels
  (
    dragForce
  );
  dragForceProps
  {
  };
\*---------------------------------------------------------------------------*/

#ifndef __DRAG_FORCE_H__
#define __DRAG_FORCE_H__

#include <functional>
#include "./force_model.h"

namespace Foam {

class dragForce : public forceModel {
 public:
  //! \brief Runtime type information
  cfdemTypeName("dragForce");

  cfdemDefineNewFunctionAdder(forceModel, dragForce);

  //! \brief Constructor
  dragForce(cfdemCloud& cloud);

  //! \brief Destructor
  ~dragForce();

  void setForce();

 private:
  //! \note subPropsDict_ should be declared in front of other members
  dictionary subPropsDict_;

  //! \brief 阻力模型名称
  std::string dragModelName_;

  //! \brief 速度场名称
  std::string velFieldName_;

  //! \brief 局部平均颗粒速度场名称
  std::string UsFieldName_;

  //! \brief 空隙率场的名称
  std::string voidFractionFieldName_;

  const volVectorField& U_;

  const volVectorField& UsField_;

  const volScalarField& voidFraction_;

  static std::hash<std::string> strHasher_;

  static const size_t DiFeliceHashValue_;

  static const size_t AbrahamHashValue_;

  static const size_t SchillerNaumannHashValue_;

  static const size_t GidaspowHashValue_;

  static const size_t SyamlalObrienHashValue_;
};

std::hash<std::string> dragForce::strHasher_;

const size_t dragForce::DiFeliceHashValue_ = strHasher_("DiFelice");

const size_t dragForce::AbrahamHashValue_ = strHasher_("Abraham");

const size_t dragForce::SchillerNaumannHashValue_ = strHasher_("Schiller-Naumann");

const size_t dragForce::GidaspowHashValue_ = strHasher_("Gidaspow");

const size_t dragForce::SyamlalObrienHashValue_ = strHasher_("Syamlal-Obrien");

}  // namespace Foam

#endif  // __DRAG_FORCE_H__
