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

#ifndef __GLOBAL_FORCE_H__
#define __GLOBAL_FORCE_H__

#include <memory>
#include <unordered_map>
#include "base/run_time_selection_tables.h"
#include "cloud/cfdem_cloud.h"
#include "sub_model/data_exchange_model/data_exchange_model.h"

namespace Foam {

//! \brief 全局力模型，封装 force model 的共享数据与共享函数
class globalForce {
 public:
  //! \brief Runtime type information
  cfdemBaseTypeName("globalForce", "");

  //! \brief Declare runtime constructor selection
  cfdemDeclareRunTimeSelection(autoPtr, globalForce, (cfdemCloud & cloud), (cloud));

  //! \brief Selector
  static autoPtr<globalForce> New(cfdemCloud& cloud, const dictionary& dict);

  //! \note base type can alse be selected.
  cfdemDefineNewFunctionAdder(globalForce, globalForce);

  globalForce(cfdemCloud& cloud);

  virtual ~globalForce();

  //! \brief 每一次耦合中，在 set force 前执行
  virtual void initBeforeSetForce();

  //! \brief 每一次耦合中，在 set force 后执行
  virtual void endAfterSetForce();

  //! \brief 更新颗粒速度
  void updatePrevParticleVelMap(const int index);

  //! \brief 更新 ddtUrHistoryMap
  void updateDDtUrHistoryMap(const int index);

  //! \brief 获取上一个耦合时间步中颗粒速度
  Foam::vector getPrevParticleVel(const int index) const;

  //! \brief 获取颗粒的历史 ddtUr
  std::vector<Foam::vector>& getDDtUrHistory(const int index);

  //! \brief 获取颗粒处背景流体速度
  virtual Foam::vector getBackgroundUfluid(const int index) const {
    FatalError << __func__ << " not implement in Foam::globalForce, please use Foam::mixGlobalForce"
               << abort(FatalError);
    return Foam::vector::zero;
  }

  //! \brief 获取颗粒处背景空隙率
  virtual double getBackgroundVoidFraction(const int index) const {
    FatalError << __func__ << " not implement in Foam::globalForce, please use Foam::mixGlobalForce"
               << abort(FatalError);
    return 1.0;
  }

  //! \brief 获取颗粒处背景流体的 ddtU
  virtual Foam::vector getBackgroundDDtU(const int index) const {
    FatalError << __func__ << " not implement in Foam::globalForce, please use Foam::mixGlobalForce"
               << abort(FatalError);
    return Foam::vector::zero;
  }

  //! \brief 获取颗粒处背景流体的涡量
  virtual Foam::vector getBackgroundVorticity(const int index) const {
    FatalError << __func__ << " not implement in Foam::globalForce, please use Foam::mixGlobalForce"
               << abort(FatalError);
    return Foam::vector::zero;
  }

  inline void resetImpParticleForce() {
    impParticleForce_ == dimensionedVector("zero", impParticleForce_.dimensions(), vector::zero);
  }

  inline void resetExpParticleForcee() {
    expParticleForce_ == dimensionedVector("zero", expParticleForce_.dimensions(), vector::zero);
  }

  inline volVectorField& impParticleForce() { return impParticleForce_; }

  inline volVectorField& expParticleForce() { return expParticleForce_; }

  inline const volVectorField& ddtU() const { return ddtU_; }

  inline const volVectorField& vorticityField() const { return vorticityField_; }

  inline const volVectorField& impParticleForce() const { return impParticleForce_; }

  inline const volVectorField& expParticleForce() const { return expParticleForce_; }

#ifdef version21
  inline const uniformDimensionedVectorField& g() const { return g_; }
#elif defined(version16ext) || defined(version15)
  inline const dimensionedVector& g() const { return g_; }
#endif

  inline const volScalarField& rhoField() const { return rho_; }

  inline const volVectorField& U() const { return U_; }

  inline const volScalarField& p() const { return p_; }

  inline const surfaceScalarField& phi() const { return phi_; }

  inline const volScalarField& voidFraction() const { return voidFraction_; }

  inline const volScalarField& volumeFraction() const { return volumeFraction_; }

 protected:
  cfdemCloud& cloud_;

  //! \note subPropsDict_ should be declared in front of other members
  dictionary subPropsDict_;

  bool verbose_;

  double GaussCoreEff_;

  //! \brief 上一个时间步中的颗粒的速度
  //! \note map: 颗粒索引 --> 颗粒速度
  std::unordered_map<int, Foam::vector> prevParticleVelMap_;

  //! \brief 记录每一个颗粒的每一个耦合时间步的 ddtUr = d(Ufluid - Up) / dt
  //! \note map: 颗粒索引 --> 历史 ddtUr
  std::unordered_map<int, std::vector<Foam::vector>> ddtUrHistoryMap_;

  volVectorField ddtU_;

  //! \brief 涡量场，vorticityField_ = fvc::curl(U_)
  volVectorField vorticityField_;

  //! \brief 颗粒隐式力的总和 [N]
  volVectorField impParticleForce_;

  //! \brief 颗粒显式力的总和 [N]
  volVectorField expParticleForce_;

  //! \brief name of the gravity field
  std::string gravityFieldName_;

  //! \brief name of the rho field
  std::string densityFieldName_;

  //! \brief 速度场名称
  std::string velFieldName_;

  //! \brief name of the finite volume pressure field
  std::string pressureFieldName_;

  //! \brief 通量场的名称
  std::string phiFieldName_;

  //! \brief 空隙率场的名称
  std::string voidFractionFieldName_;

  //! \brief name of the finite volume voidfraction field
  std::string volumeFractionFieldName_;

#ifdef version21
  const uniformDimensionedVectorField& g_;
#elif defined(version16ext) || defined(version15)
  const dimensionedVector& g_;
#endif

  const volScalarField& rho_;

  const volVectorField& U_;

  const volScalarField& p_;

  const surfaceScalarField& phi_;

  const volScalarField& voidFraction_;

  const volScalarField& volumeFraction_;
};

}  // namespace Foam

#endif  // __GLOBAL_FORCE_H__
