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
  forceSubModel
\*---------------------------------------------------------------------------*/

#ifndef __FORCE_SUB_MODEL_H__
#define __FORCE_SUB_MODEL_H__

#include "cloud/cfdem_cloud.h"

namespace Foam {

//! \brief force switch type enum
enum ESwitch {
  /*!
   * \brief kTreatForceExplicitInMomEquation
   *   true - 在CFD的动量方程中，耦合力都为显式力
   *   false - 在CFD的动量方程中，耦合力都为隐式力 (default)
   */
  kTreatForceExplicitInMomEquation = 0,
  /*!
   * \brief kTreatForceBothCFDAndDEM
   *   true - 即 CFD 和 DEM 都考虑耦合力 (default)
   *   false - 仅在 DEM 中考虑耦合力，即仅考虑流体对颗粒的作用，不考虑颗粒对流体的反作用
   *           在 resolved 方法中，比如 cfdemSolverIB，ArchimedesIB force and ShirgaonkarIB
   *           会直接设置 treatForceDEM 为 true，因为颗粒对流体的作用力是通过虚拟域方法得到
   */
  kTreatForceBothCFDAndDEM,
  /*!
   * \brief kTreatDEMForceImplicit
   *   true - 在每个耦合时间步，流体的速度和阻力系数都被传递到 DEM 中，从而在每个 DEM 时间步中，
   *          使用上一个耦合时间步中的阻力系数和流体速度，与当前颗粒速度一起计算颗粒受到的阻力
   *   false - 在每个耦合时间步中，流体对颗粒的阻力被传递到 DEM 中，并且在接下来
   *           的 DEM 时间步中，这个力保持不变，直到下个耦合时间步 (default)
   */
  kTreatDEMForceImplicit,
  /*!
   * \brief kVerbose
   *   true - 调试信息输出到屏幕
   *   false - 调试信息不输出到屏幕 (default)
   */
  kVerbose,
  /*!
   * \brief kInterpolation
   *   true - 将欧拉场使用插值模型插值到拉格朗日场
   *   false - 只使用网格中心的值 (default)
   */
  kInterpolation,
  /*!
   * \brief kScalarViscosity
   *   true - 使用用户自定义动力粘度 nu 进行阻力计算，CFD 计算仍然使用 transportDict 指定的动力粘度
   *   false - 不使用用户自定义动力粘度 nu 进行阻力计算 (default)
   */
  kScalarViscosity
};

//! \brief force type enum
enum EForceType {
  //! \brief force used for unresolved method
  kUnResolved = 0,
  //! \brief force used for semi-resolved method
  kSemiResolved,
  //! \brief force used for resolved method
  kResolved,
  //! \brief force used for mix method
  kMix
};

/*!
 * \brief The force sub model is designed to hold the settings a force model can have.
 *        For now it handles the kTreatForceExplicitInMomEquation, kTreatForceBothCFDAndDEM,
 *        kTreatDEMForceImplicit and kVerbose option.
 */
class forceSubModel {
 public:
  //! \brief Constructor
  forceSubModel(cfdemCloud& cloud, forceModel& forceModel, const dictionary& subPropsDict);

  //! \brief Destructor
  ~forceSubModel();

  class Switches {
   public:
    static const int kNum;
    static const char* kNameList[];
    Switches() : value_(0) {}
    inline void setTrue(ESwitch value) {
      switch (value) {
        case kTreatForceExplicitInMomEquation:
          value_ |= (1 << kTreatForceExplicitInMomEquation);
          break;
        case kTreatForceBothCFDAndDEM:
          value_ |= (1 << kTreatForceBothCFDAndDEM);
          break;
        case kTreatDEMForceImplicit:
          value_ |= (1 << kTreatDEMForceImplicit);
          break;
        case kVerbose:
          value_ |= (1 << kVerbose);
          break;
        case kInterpolation:
          value_ |= (1 << kInterpolation);
          break;
        case kScalarViscosity:
          value_ |= (1 << kScalarViscosity);
          break;
        default:
          FatalError << "Error: illegal switch enum: " << value << " in forceSubModel" << abort(FatalError);
      }
    }
    inline bool isTrue(ESwitch value) const {
      switch (value) {
        case kTreatForceExplicitInMomEquation:
          return value_ & (1 << kTreatForceExplicitInMomEquation);
        case kTreatForceBothCFDAndDEM:
          return value_ & (1 << kTreatForceBothCFDAndDEM);
        case kTreatDEMForceImplicit:
          return value_ & (1 << kTreatDEMForceImplicit);
        case kVerbose:
          return value_ & (1 << kVerbose);
        case kInterpolation:
          return value_ & (1 << kInterpolation);
        case kScalarViscosity:
          return value_ & (1 << kScalarViscosity);
        default:
          FatalError << "Error: illegal switch enum: " << value << " in forceSubModel" << abort(FatalError);
      }
      return false;
    }
    inline bool isFalse(ESwitch value) const { return !isTrue(value); }

   private:
    unsigned int value_;
  };

  /*!
   * \param index                  <[in] 颗粒索引
   * \param dragTot                <[in] 索引为 index 的颗粒受到的总阻力
   * \param dragEx                 <[in] 索引为 index 的颗粒受到的显式阻力
   * \param Ufluid = vector::zero  <[in] 索引为 index 的颗粒中心处流体速度(可以指定是否使用插值模型计算)
   * \param scalar Cd = 0          <[in] 颗粒阻力系数
   */
  void partToArray(const int& index, const Foam::vector& dragTot, const Foam::vector& dragEx,
                   const Foam::vector& Ufluid = Foam::vector::zero, scalar Cd = scalar(0)) const;

  /*!
   * \param index                  <[in] 颗粒索引
   * \param torque = vector::zero  <[in] 索引为 index 的颗粒受到的力矩
   */
  void addTorque(int index, const Foam::vector& torque = Foam::vector::zero) const;

  //! \brief read switches from force model dictionary
  void readSwitches();

  /*!
   * \brief check switches for different type of force
   * \param forceType force type, Eg: kUnResolved
   */
  void checkSwitches(EForceType forceType) const;

  /*!
   * \brief 计算 IB drag，用于计算 ShirgaonkarIBModel and mixShirgaonkarIBModel
   *        dimensionSet(1, -2, -2, 0, 0)
   * \param U 速度场
   * \param p 压力场
   * \return IB drag
   */
  volVectorField IBDrag(const volVectorField& U, const volScalarField& p) const;

  inline const volScalarField& rhoField() const { return rho_; }

  inline const volScalarField& nuField() const {
#ifdef compre
    nu_ = cloud_.turbulence().mu() / rho_;
    return nu_;
#else
    if (switches_.isTrue(kScalarViscosity)) {
      return nu_;
    } else {
      return cloud_.turbulence().nu();
    }
#endif
  }

  inline const volScalarField& muField() const {
#ifdef compre
    return cloud_.turbulence().mu();
#else
    if (switches_.isTrue(kScalarViscosity)) {
      return nu_ * rho_;
    } else {
      return cloud_.turbulence().nu() * rho_;
    }
#endif
  }

  inline bool treatForceExplicitInMomEquation() const { return switches_.isTrue(kTreatForceExplicitInMomEquation); }

  inline bool treatForceBothCFDAndDEM() const { return switches_.isTrue(kTreatForceBothCFDAndDEM); }

  inline bool treatDEMForceImplicit() const { return switches_.isTrue(kTreatDEMForceImplicit); }

  inline bool verbose() const { return switches_.isTrue(kVerbose); }

  inline bool interpolation() const { return switches_.isTrue(kInterpolation); }

  inline bool scalarViscosity() const { return switches_.isTrue(kScalarViscosity); }

 protected:
  cfdemCloud& cloud_;

  forceModel& forceModel_;

  const dictionary subPropsDict_;

  Switches switches_;

  const std::string densityFieldName_;

  //! \brief 密度
  const volScalarField& rho_;

  //! \brief 自定义动力粘度场
  volScalarField nu_;
};

}  // namespace Foam

#endif  // __FORCE_SUB_MODEL_H__
