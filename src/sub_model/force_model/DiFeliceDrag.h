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
  interaction forces) acting on each DEM particle. The DiFeliceDrag model
  is a model that calculates the particle based drag force following the
  correlation of Di Felice (see Zhou et al. (2010), JFM).

Syntax
  forceModels
  (
    DiFeliceDrag
  );
  DiFeliceDragProps
  {
    velFieldName "U";
    voidFractionFieldName "voidFraction";
    granVelFieldName "Us";
    interpolation switch1;
    voidFractionInterpolationType "type1";
    UInterpolationType "type2";
    suppressProbe       switch2;
    scale               scalar1;
    scaleDrag           scalar2;
    treatForceExplicit  switch3;
    implForceDEM        switch4;
    verbose             switch5;
    scalarViscosity     switch6;
    nu                  scalar3;
  };
\*---------------------------------------------------------------------------*/

#ifndef __DIFELICE_DRAG_H__
#define __DIFELICE_DRAG_H__

#include "./force_model.h"

namespace Foam {

class DiFeliceDrag : public forceModel {
 public:
  //! \brief Runtime type information
  cfdemTypeName("DiFeliceDrag");

  cfdemDefineNewFunctionAdder(forceModel, DiFeliceDrag);

  //! \brief Constructor
  DiFeliceDrag(cfdemCloud& cloud);

  //! \brief Destructor
  ~DiFeliceDrag();

  void setForce();

 private:
  //! \note subPropsDict_ should be declared in front of other members
  dictionary subPropsDict_;

  //! \brief 速度场名称
  std::string velFieldName_;

  //! \brief 局部平均颗粒速度场名称
  std::string UsFieldName_;

  //! \brief 空隙率场的名称
  std::string voidFractionFieldName_;

  const volVectorField& U_;

  const volVectorField& UsField_;

  const volScalarField& voidFraction_;
};

}  // namespace Foam

#endif  // __DIFELICE_DRAG_H__
