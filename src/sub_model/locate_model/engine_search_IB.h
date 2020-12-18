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
  The locateModel “engineIB” locates the CFD cell and cellID corresponding
  to a given position. This locate model is especially designed for parallel
  immersed boundary method. Each particle is represented by “satellite points”
  if it is distributed over several processors.

  The engineSearchIB locate Model can be used with different settings to
  use different algorithms:
  - treeSearch false; will execute some geometric (linear) search using the
    last known cellID (recommended)
  - treeSearch true; will use a recursive tree structure to find the cell.

  This model is a modification of the engine search model. Instead of using
  the centre-cell as starting point for the engine search, further satellite
  points located on the surface of the sphere are checked. This ensures that
  (parts of) spheres can be located even when their centre is on another
  processor. This is especially important for parallel computations, when a
  sphere is about to move from one processor to another.

Syntax
  locateModel engineIB;
  engineIBProps
  {
    treeSearch switch1;
    zSplit value1;
    xySplit value2;
  }

Restrictions
  Only for immersed boundary solvers.

Class
  engineSearchIB
\*---------------------------------------------------------------------------*/

#include "./engine_search.h"

namespace Foam {

class engineSearchIB : public engineSearch {
 public:
  //! \brief Runtime type information
  cfdemTypeName("engineIB");

  cfdemDefineNewFunctionAdder(locateModel, engineSearchIB);

  //! \brief Constructor
  engineSearchIB(cfdemCloud& cloud, const std::string& typeName = cTypeName());

  //! \brief Destructor
  ~engineSearchIB();

  /*!
   * \brief use search engine to get cell id of particle center
   * \param findCellIDs 颗粒覆盖网格的编号
   */
  void findCell(const base::CITensor1& findCellIDs) const;

 private:
  /*!
   * \brief 判断 pos 是否位于长方体区域中
   * \param pos 颗粒中心位置
   * \param offsetValue 位置偏移量
   */
  bool isInsideRectangularDomain(const Foam::vector& pos, double offsetValue) const;

  /*!
   * \brief generate satellite point according to index
   * \param index satellite point index
   */
  Foam::vector generateSatellitePoint(int index) const;

  /*!
   * \brief get satellite point position
   * \param index particle index
   * \param satellitePointIndex satellite point index
   */
  inline Foam::vector getSatellitePointPos(int index, int satellitePointIndex) const;

  //! \brief 如果网格更新，则调用该函数修正 searchEngine_ 以及重新设置 boundBox
  virtual void correctSearchEngine() {
    engineSearch::correctSearchEngine();
    // searchEngine_.correct();
    boundBoxPtr_.reset(new boundBox(cloud_.mesh().points(), false));
    if (verbose_) {
      Pout << "Min bounds (x, y, z): " << boundBoxPtr_().min() << endl;
      Pout << "Max bounds (x, y, z): " << boundBoxPtr_().max() << endl;
    }
  }

 private:
  bool verbose_;

  int zSplit_;

  int xySplit_;

  int numberOfSatellitePoints_;

  double coef_;

  autoPtr<boundBox> boundBoxPtr_;

  std::vector<Foam::vector> satellitePoints_;
};

}  // namespace Foam
