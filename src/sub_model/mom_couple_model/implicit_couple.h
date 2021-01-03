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
  The implicitCouple-model is a momCoupleModel model providing an implicit
  momentum source term for the CFD solver.

Syntax
  Defined in couplingProperties dictionary:
  momCoupleModels
  (
    implicitCouple
  );
  implicitCoupleProps
  {
    velFieldName "U";
    granVelFieldName "Us";
    voidFractionFieldName "voidFraction";
    minAlphaP number;
  }
  U = name of the finite volume fluid velocity field
  Us = name of the finite volume granular velocity field
  voidFraction = name of the finite volume voidFraction field
  number = minimum value for local particle volume fraction to calculate
           the exchange filed (default Foam::SMALL)

Restrictions
  Only for solvers that include implicit momentum exchange.

Class
  Foam::implicitCouple
\*---------------------------------------------------------------------------*/

#ifndef __IMPLICIT_COUPLE_H__
#define __IMPLICIT_COUPLE_H__

#include "./mom_couple_model.h"

namespace Foam {

class implicitCouple : public momCoupleModel {
 public:
  //! \brief Runtime type information
  cfdemTypeName("implicitCouple");

  cfdemDefineNewFunctionAdder(momCoupleModel, implicitCouple);

  //! \brief Constructor
  implicitCouple(cfdemCloud& cloud);

  //! \brief Destructor
  ~implicitCouple();

  //! \brief reset mom field
  void resetMomSourceField();

  //! \brief get implicit momentum source field
  tmp<volScalarField> impMomSource();

 private:
  //! \note subPropsDict_ should be declared in front of other members
  dictionary subPropsDict_;

  //! \brief 流体速度场的名称
  std::string velFieldName_;

  //! \brief 局部平均颗粒速度场的名称
  std::string granVelFieldName_;

  //! \brief 空隙率场的名称
  std::string voidFractionFieldName_;

  //! \brief 速度场的常引用
  const volVectorField& U_;

  //! \brief 局部平均颗粒速度场的常引用
  const volVectorField& Us_;

  //! \brief 空隙率场的常引用
  const volScalarField& alpha_;

  //! @brief 动量交换场的最大值
  double KslLimit_;

  /*!
   * \brief 动量交换场
   * \note Ksl = 单位体积上颗粒对流体的作用力 / (Us - U)，单位为 (N/m^3) / (m/s) = kg / (m^3 * s)
   *             dimensionedScalar("zero", dimensionSet(1,-3,-1,0,0), 0)
   */
  volScalarField KslPrev_;

  volScalarField KslNext_;
};

}  // namespace Foam

#endif  // __IMPLICIT_COUPLE_H__
