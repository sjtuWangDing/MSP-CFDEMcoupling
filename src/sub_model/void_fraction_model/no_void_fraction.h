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
  The noVoidFraction voidFraction model is a dummy model and has no
  physical meaning.

Syntax
  voidfractionModel noVoidFraction;

Class
  noVoidFraction
\*---------------------------------------------------------------------------*/

#ifndef __NO_VOID_FRACTION_H__
#define __NO_VOID_FRACTION_H__

#include "./void_fraction_model.h"

namespace Foam {

class noVoidFraction : public voidFractionModel {
 public:
  cfdemTypeName("noVoidFraction");

  cfdemDefineNewFunctionAdder(voidFractionModel, noVoidFraction);

  //! \brief Constructor
  noVoidFraction(cfdemCloud& cloud);

  //! \brief Destructor
  virtual ~noVoidFraction();
};

}  // namespace Foam

#endif  // __NO_VOID_FRACTION_H__
