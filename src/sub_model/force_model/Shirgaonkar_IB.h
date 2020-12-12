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
  interaction forces) acting on each DEM particle. The ShirgaonkarIB model
  calculates the drag force (viscous and pressure force) acting on each
  particle in a resolved manner (see Shirgaonkar et al. (2009): “A new
  mathematical formulation and fast algorithm for fully resolved simulation
  of self-propulsion”, Journal of Comp. Physics). This model is only suited
  for resolved CFD-DEM simulations where the particle is represented by
  immersed boundary method.

Syntax
  forceModels
  (
    ShirgaonkarIB
  );
  ShirgaonkarIBProps
  {
    velFieldName "U";
    pressureFieldName "p";
    verbose true;
    treatForceExplicit true;
  };

Restrictions
  Only for immersed boundary solvers.
\*---------------------------------------------------------------------------*/

#ifndef __SHIRGAONKAR_IB_H__
#define __SHIRGAONKAR_IB_H__

#include "./force_model.h"

namespace Foam {

class ShirgaonkarIB : public forceModel {
 public:
  //! \brief Runtime type information
  cfdemTypeName("ShirgaonkarIB")

      cfdemDefineNewFunctionAdder(forceModel, ShirgaonkarIB)

      //! \brief Constructor
      ShirgaonkarIB(cfdemCloud& cloud);

  //! \brief Destructor
  ~ShirgaonkarIB();

  void setForce();

 private:
  //! \note subPropsDict_ should be declared in front of other members
  dictionary subPropsDict_;

  //! \brief name of the finite volume fluid velocity field
  std::string velFieldName_;

  //! \brief name of the finite volume pressure field
  std::string pressureFieldName_;

  const volVectorField& U_;

  const volScalarField& p_;

  bool useTorque_;
};

}  // namespace Foam

#endif  // __SHIRGAONKAR_IB_H__
