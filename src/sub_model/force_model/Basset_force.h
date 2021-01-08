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

#ifndef __BASSET_FORCE_H__
#define __BASSET_FORCE_H__

#include <unordered_map>
#include "./force_model.h"

namespace Foam {

class BassetForce : public forceModel {
 public:
  //! \brief Runtime type information
  cfdemTypeName("BassetForce");

  cfdemDefineNewFunctionAdder(forceModel, BassetForce);

  //! \brief Constructor
  BassetForce(cfdemCloud& cloud);

  //! \brief Destructor
  ~BassetForce();

  void setForce();

 protected:
  void setMiddleParticleForceKernel(Foam::vector& BstForce, const int index, const volScalarField& rhoField,
                                    const volScalarField& nuField, const volVectorField& ddtUField);

  void updatePrevUp(const int index);

  void updateDDtUrHistory(const int index);

 private:
  //! \note subPropsDict_ should be declared in front of other members
  dictionary subPropsDict_;

  Foam::vector initUr_;

  //! \brief 记录每一个颗粒的上一个耦合时间步的速度
  std::unordered_map<int, Foam::vector> prevUpMap_;

  //! \brief 记录每一个颗粒的历史 (dU / dt - dV / dt)
  std::unordered_map<int, std::vector<Foam::vector>> DDtUrHistoryMap_;

  //! \brief name of the finite volume fluid velocity field
  std::string velFieldName_;

  //! \brief 空隙率场的名称
  std::string voidFractionFieldName_;

  //! \brief name of the phi field
  std::string phiFieldName_;

  const volVectorField& U_;

  const volScalarField& voidFraction_;

  const surfaceScalarField& phi_;
};

}  // namespace Foam

#endif  // __BASSET_FORCE_H__
