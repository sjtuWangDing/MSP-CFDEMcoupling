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

#include "./dense.h"

namespace Foam {

cfdemDefineTypeName(dense);

cfdemCreateNewFunctionAdder(averagingModel, dense);

//! \brief Constructor
dense::dense(cfdemCloud& cloud) : averagingModel(cloud) {}

//! \brief Destructor
dense::~dense() {}

/*!
 * \brief 计算局部平均矢量场
 * \param valueField   <[in, out] 需要被局部平均化场
 * \param weightField  <[in, out] 权重系数平均化场
 * \param value        <[in] 用于局部平均化的颗粒数据(lagrange value)
 * \param weight       <[in] 用于局部平均化的权重系数(lagrange value)
 */
void dense::setVectorFieldAverage(volVectorField& valueField, volScalarField& weightField,
                                  const base::CDExTensor2& value, const std::vector<base::CDTensor1>& weight) {
  CHECK_EQ(value.size(1), 3) << __func__ << ": vector field's dimension must equal to 3";
  for (int index = 0; index < cloud_.numberOfParticles(); index++) {
    if (cloud_.checkFineParticle(index) || cloud_.checkMiddleParticle(index)) {
      // get value vector
      Foam::vector valueVec(value[index][0], value[index][1], value[index][2]);
      for (int subCell = 0; subCell < cloud_.particleOverMeshNumber()[index]; ++subCell) {
        int cellID = cloud_.cellIDs()[index][subCell];
        if (cellID >= 0) {  // cell Found
          double weightP = weight[index][subCell];
          if (fabs(weightField[cellID]) < Foam::SMALL) {
            // first entry in this cell
            valueField[cellID] = valueVec;
            weightField[cellID] = weightP;
          } else {
            // not first entry in this cell
            valueField[cellID] =
                (valueField[cellID] * weightField[cellID] + valueVec * weightP) / (weightField[cellID] + weightP);
            weightField[cellID] += weightP;
          }
        }
      }
    }
  }
  valueField.correctBoundaryConditions();
  base::MPI_Info("dense: setVectorFieldAverage - done", true);
}

}  // namespace Foam
