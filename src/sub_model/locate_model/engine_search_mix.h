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

#ifndef __ENGINE_SEARCH_MIX_H__
#define __ENGINE_SEARCH_MIX_H__

#include "./engine_search_IB.h"

namespace Foam {

class engineSearchMix : public engineSearchIB {
 public:
  //! \brief Runtime type information
  cfdemTypeName("engineMix");

  cfdemDefineNewFunctionAdder(locateModel, engineSearchMix);

  //! \brief Constructor
  engineSearchMix(cfdemCloud& cloud, const std::string& typeName = cTypeName());

  //! \brief Destructor
  ~engineSearchMix();

  /*!
   * \brief use search engine to get cell id of particle center
   * \param findCellIDs 颗粒中心覆盖网格的编号
   */
  void findCell(const base::CITensor1& findCellIDs) const;

  /*!
   * \brief use search engine to get id of cell which covered by processor
   * \param findMpiCellIDs 颗粒覆盖当前求解器上任意一个网格的编号
   * \param scale 颗粒半径尺度系数
   */
  void findMpiCell(const base::CITensor1& findMpiCellIDs) const;

  /*!
   * \brief use search engine to get id of cell which covered by processor
   * \param findExpandedCellIDs 被扩展的颗粒覆盖当前求解器上任意一个网格的编号
   */
  void findExpandedCell(const base::CITensor1& findExpandedCellIDs, const double scale) const;

 protected:
  //! \brief 每个颗粒中心只可能位于一个网格中，但是如果颗粒位于处理器边界上，则颗粒会被边界两边的处理器都捕获到
  //!   所以该函数会确保所有颗粒只被一个处理器捕获到
  void uniqueFindCell(const base::CITensor1& findCellIDs) const;

  void findMpiCellKernel(const int index, const base::CITensor1& findMpiCellIDs, const double scale = 1.0) const;
};

}  // namespace Foam

#endif  // __ENGINE_SEARCH_MIX_H__