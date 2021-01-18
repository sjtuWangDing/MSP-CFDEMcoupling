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
  The locateModel “engine” locates the CFD cell and cellID corresponding
  to a given position.

Syntax
  locateModel engine;
  engineProps
  {
    treeSearch switch1;
  }

Class
  engineSearch
\*---------------------------------------------------------------------------*/

#include "./locate_model.h"
#include "meshSearch.H"

namespace Foam {

class engineSearch : public locateModel {
 public:
  //! \brief Runtime type information
  cfdemTypeName("engine");

  cfdemDefineNewFunctionAdder(locateModel, engineSearch);

  //! \brief Constructor
  engineSearch(cfdemCloud& cloud, const std::string& derivedTypeName = cTypeName());

  //! \brief Destructor
  ~engineSearch();

  /*!
   * \brief use search engine to get cell id of particle center
   * \param findCellIDs 颗粒覆盖网格的编号
   */
  void findCell(const base::CITensor1& findCellIDs) const;

  //! \brief 每个颗粒中心只可能位于一个网格中，但是如果颗粒位于处理器边界上，则颗粒会被边界两边的处理器都捕获到
  //!   所以该函数会确保所有颗粒只被一个处理器捕获到
  void uniqueFindCell(const base::CITensor1& findCellIDs) const;

  /*!
   * \brief use search engine to get id of cell which covered by processor
   * \param findMpiCellIDs 颗粒覆盖当前求解器上任意一个网格的编号
   */
  void findMpiCell(const base::CITensor1& findMpiCellIDs) const {
    FatalError << __func__ << " not implemented in engineSearch model, please use enginehMix model\n"
               << abort(FatalError);
  }

  /*!
   * \brief use search engine to get id of cell which covered by processor
   * \param findExpandedCellIDs 被扩展的颗粒覆盖当前求解器上任意一个网格的编号
   */
  void findExpandedCell(const base::CITensor1& findExpandedCellIDs, const double scale) const {
    FatalError << __func__ << " not implemented engineSearch model, please use enginehMix model\n" << abort(FatalError);
  }

  /*!
   * \brief use search engine to get cell id of certain vector
   * \param position 颗粒中心的位置
   * \param oldCellID old cell ID
   */
  label findSingleCell(const Foam::vector& position, label oldCellID) const;

  //! \brief 如果网格更新，则调用该函数修正 searchEngine_
  virtual void correctSearchEngine() { searchEngine_.correct(); }

  inline Foam::meshSearch& searchEngine() { return searchEngine_; }

 protected:
  dictionary subPropsDict_;

  //! \brief 是否开启 octree 搜索，default - true
  bool treeSearch_;

  Foam::meshSearch searchEngine_;
};

}  // namespace Foam
