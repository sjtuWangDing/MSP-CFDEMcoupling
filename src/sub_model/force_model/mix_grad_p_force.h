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

#ifndef __MIX_GRAD_P_FORCE_H__
#define __MIX_GRAD_P_FORCE_H__

#include "./force_model.h"

namespace Foam {

class mixGradPForce : public forceModel {
 public:
  //! \brief Runtime type information
  cfdemTypeName("mixGradPForce");

  cfdemDefineNewFunctionAdder(forceModel, mixGradPForce);

  //! \brief Constructor
  mixGradPForce(cfdemCloud& cloud);

  //! \brief Destructor
  ~mixGradPForce();

  void setForce();

 protected:
  void setForceKernel(const int index, Foam::vector& mixGradPF, std::once_flag& onceOp);

  Foam::vector getBackgroundGradP(const int index, const int findCellID) const;

 private:
  //! \note subPropsDict_ should be declared in front of other members
  dictionary subPropsDict_;

  const volScalarField& p_;

  const volVectorField& gradPField_;
};

}  // namespace Foam

#endif  // __MIX_GRAD_P_FORCE_H__
