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
  This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
  and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

Class
  Foam::twoWayMPI
\*---------------------------------------------------------------------------*/

#ifndef __TWO_WAY_MPI_H__
#define __TWO_WAY_MPI_H__

#include "mpi.h"
#include <lammps.h>
#include <input.h>
#include <atom.h>
#include <library.h>
#include <error.h>
#include <update.h>
#include <library_cfd_coupling.h>

#include "./data_exchange_model.h"

namespace Foam {

class twoWayMPI : public dataExchangeModel {

public:

  //! \brief Runtime type information
  cfdemTypeName("twoWayMPI")

  cfdemDefineNewFunctionAdder(dataExchangeModel, twoWayMPI)

  //! \brief Constructor
  twoWayMPI(cfdemCloud& cloud);

  //! \brief Destructor
  ~twoWayMPI();

  //! \return 当前耦合时间步中颗粒的数量
  int couple();

  // void destroy(double**& array);
  // void destroy(double**& array, int);
  // void allocateArray(double**& array, double initVal, int width, int length);

  // //! \brief Allocate and destroy for 2-D int array
  // void destroy(int**& array);
  // void destroy(int**& array, int);
  // void allocateArray(int**& array, int initVal, int width, int length);
  // void allocateArray(int**& array, int initVal, int width, const char* length  = "nparticles");

  // //! \brief Allocate and destroy for 1-D double array
  // void destroy(double*& array);
  // void allocateArray(double*& array, double initVal, int length);

  // //! \brief Allocate and destroy for 1-D int array
  // void destroy(int*& array);
  // void allocateArray(int*& array, int initVal, int length);

  /*!
   * \brief 从 LAMMPS 中获取数据（数据类型为 double）
   *        Eg: getData("x", "vector-atom", array);
   * \param dataName 数据名称
   * \param dataType 数据类型
   * \param array 内存地址（LIGGGHTS code 要求 array 为二级指针）
   */
  void getData(const std::string& dataName,
               const std::string& dataType,
               double** array) const {
    // 虽然这里将 array 转为了 void*，但是在 LIGGGHTS code 中会将其转为 double**
    void* p = static_cast<void*>(array);
    data_liggghts_to_of(dataName.c_str(), dataType.c_str(), lmp_, p, (const char*)"double");
  }

  /*!
   * \brief 从 LAMMPS 中获取数据（数据类型为 int）
   * \param dataName 数据名称
   * \param dataType 数据类型
   * \param array 内存地址（LIGGGHTS code 要求 array 为二级指针）
   */
  void getData(const std::string& dataName,
               const std::string& dataType,
               int** array) const {
    // 虽然这里将 array 转为了 void*，但是在 LIGGGHTS code 中会将其转为 int**
    void* p = static_cast<void*>(array);
    data_liggghts_to_of(dataName.c_str(), dataType.c_str(), lmp_, p, (const char*)"int");
  }

  /*!
   * \brief 传递数据到 LAMMPS 中（数据类型为 double）
   *        Eg: giveData("dragforce", "vector-atom", array);
   * \param dataName 数据名称
   * \param dataType 数据类型
   * \param array 内存地址（LIGGGHTS code 要求 array 为二级指针）
   */
  void giveData(const std::string& dataName,
                const std::string& dataType,
                double** field) const {
    data_of_to_liggghts(dataName.c_str(), dataType.c_str(), lmp_, (void*)field, (const char*)"double");
  }

  /*!
   * \brief 传递数据到 LAMMPS 中（数据类型为 int）
   * \param dataName 数据名称
   * \param dataType 数据类型
   * \param array 内存地址（LIGGGHTS code 要求 array 为二级指针）
   */
  void giveData(const std::string& dataName,
                const std::string& dataType,
                int** field) const {
    data_of_to_liggghts(dataName.c_str(), dataType.c_str(), lmp_, (void*)field, (const char*)"int");
  }

  //! \brief Allocate for 2-D double array using liggghts interface
  void liggghtsAllocate(double**& array, int length, int width, double initVal = 0.0) const {
    allocate_external_double(array, width, length, initVal, lmp_);
  }

  //! \brief Allocate 2-D int array using liggghts interface
  void liggghtsAllocate(int**& array, int length, int width, int initVal = 0) const {
    allocate_external_int(array, width, length, initVal, lmp_);
  }

private:

  LAMMPS_NS::LAMMPS* lmp_;

  MPI_Comm lgs_comm_;

  dictionary subPropsDict_;
};

} // namespace Foam

#endif // __TWO_WAY_MPI_H__
