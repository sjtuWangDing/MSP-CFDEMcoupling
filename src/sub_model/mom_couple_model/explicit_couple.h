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
  The explicitCouple-model is a momCoupleModel model providing an explicit
  momentum source term for the CFD solver.
\*---------------------------------------------------------------------------*/

#ifndef __EXPLICIT_COUPLE_H__
#define __EXPLICIT_COUPLE_H__

#include "./mom_couple_model.h"

namespace Foam {

class explicitCouple : public momCoupleModel {
 public:
  //! \brief Runtime type information
  cfdemTypeName("explicitCouple");

  cfdemDefineNewFunctionAdder(momCoupleModel, explicitCouple);

  //! \brief Constructor
  explicitCouple(cfdemCloud& cloud);

  //! \brief Destructor
  ~explicitCouple();

  //! \brief reset mom field
  void resetMomSourceField();

  //! \brief get explicit momentum source field
  tmp<volVectorField> expMomSource();

 private:
  //! \note subPropsDict_ should be declared in front of other members
  dictionary subPropsDict_;

  //! @brief 显式力场的最大值
  Foam::vector expForceLimit_;

  /*!
   * \brief 显式力场
   * \note expForce = 单位体积上颗粒对流体的作用力，单位为 kg / (m^2 * s^2)
   */
  volVectorField expForcePrev_;

  volVectorField expForceNext_;
};

}  // namespace Foam

#endif  // __EXPLICIT_COUPLE_H__
