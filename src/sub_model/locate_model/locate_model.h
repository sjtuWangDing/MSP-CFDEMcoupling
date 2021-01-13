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
  The locateModel is the base class for models which search for the CFD
  cell and cellID corresponding to a position. In general it is used to
  find the cell a particle is located in.

Syntax
  locateModel modelName;

Class
  locateModel
\*---------------------------------------------------------------------------*/

#ifndef __LOCAL_MODEL_H__
#define __LOCAL_MODEL_H__

#include "base/run_time_selection_tables.h"
#include "cloud/cfdem_cloud.h"

namespace Foam {

//! \brief locateModel
class locateModel {
 public:
  //! \brief Runtime type information
  cfdemBaseTypeName("locateModel", "");

  //! \brief Declare runtime constructor selection
  cfdemDeclareRunTimeSelection(autoPtr, locateModel, (cfdemCloud & cloud), (cloud));

  //! \brief Selector
  static autoPtr<locateModel> New(cfdemCloud& cloud, const dictionary& dict);

  //! \brief Constructor
  locateModel(cfdemCloud& cloud);

  //! \brief Destructor
  virtual ~locateModel();

  /*!
   * \brief use search engine to get cell id of particle center
   * \param findCellIDs 颗粒中心覆盖网格的编号
   */
  virtual void findCell(const base::CITensor1& findCellIDs) const = 0;

  virtual label findSingleCell(const Foam::vector& position, label oldCellID) const = 0;

  virtual void correctSearchEngine() {
    FatalError << "locateModel::correctSearchEngine(): in current locate model has no search engine."
               << abort(FatalError);
  }

 protected:
  cfdemCloud& cloud_;
};

}  // namespace Foam

#endif  // __LOCAL_MODEL_H__
