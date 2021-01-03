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
  momentum source term for the CFD solver and additionally it superposes an
  additional source field which can be set via the function setSourceField.

Syntax
  Defined in couplingProperties dictionary:
  momCoupleModels
  (
    explicitCouple
  );
  explicitCoupleProps
  {
    fLimit vector;
  }
  vector = limiter vector for explicit force term (default (1e10,1e10,1e10) )

Restrictions
  Only for solvers that include explicit momentum exchange.

Class
  Foam::explicitCouple
\*---------------------------------------------------------------------------*/

#ifndef __EXPLICIT_COUPLE_H__
#define __EXPLICIT_COUPLE_H__

#include "./mom_couple_model.h"

namespace Foam {

class explicitCouple : public momCoupleModel {};

}  // namespace Foam

#endif  // __EXPLICIT_COUPLE_H__
