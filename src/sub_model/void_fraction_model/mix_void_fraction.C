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

#include "./IB_void_fraction.h"
#include "./divided_void_fraction.h"
#include "./mix_void_fraction.h"
#include "sub_model/force_model/global_force.h"
#include "sub_model/locate_model/locate_model.h"

namespace Foam {

cfdemDefineTypeName(mixVoidFraction);

cfdemCreateNewFunctionAdder(voidFractionModel, mixVoidFraction);

//! \brief Constructor
mixVoidFraction::mixVoidFraction(cfdemCloud& cloud)
    : voidFractionModel(cloud),
      subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")),
      verbose_(false),
      useGuassVoidFractionForMiddleParticle_(false),
      GaussKernelBandWidth_(0.0),
      GaussKernelScale_(0.0),
      alphaMin_(0.0),
      tooMuch_(0.0) {
  weight_ = subPropsDict_.lookupOrDefault<double>("weight", 1.0);
  porosity_ = subPropsDict_.lookupOrDefault<double>("porosity", 1.0);
  verbose_ = subPropsDict_.lookupOrDefault<bool>("verbose", false);
  useGuassVoidFractionForMiddleParticle_ =
      subPropsDict_.lookupOrDefault<bool>("useGuassVoidFractionForMiddleParticle", false);
  if (useGuassVoidFractionForMiddleParticle_) {
    GaussKernelBandWidth_ = subPropsDict_.lookupOrDefault<double>("GaussKernelBandWidth", 0.0);
    GaussKernelScale_ = subPropsDict_.lookupOrDefault<double>("GaussKernelScale", 2.0);
    if (GaussKernelBandWidth_ < Foam::SMALL) {
      FatalError << __func__ << ": GaussKernelBandWidth shloud be >= 2.0 * dp(diamter of middle particle)."
                 << abort(FatalError);
    }
    if (GaussKernelScale_ < 2.0 - Foam::SMALL) {
      FatalError << __func__ << ": GaussKernelScale shloud be >= 2.0." << abort(FatalError);
    }
  }
  alphaMin_ = subPropsDict_.lookupOrDefault<double>("alphaMin", 0.0);
  // 单个颗粒覆盖最多网格数量
  maxCellsNumPerCoarseParticle_ = subPropsDict_.lookupOrDefault<int>("maxCellsNumPerCoarseParticle", 1000);
  if (alphaMin_ > 1 || alphaMin_ < 0.01) {
    FatalError << __func__ << ": alphaMin shloud be > 1 and < 0.01." << abort(FatalError);
  }
}

//! \brief Destructor
mixVoidFraction::~mixVoidFraction() {}

//! \brief 输出空隙率相关信息
void mixVoidFraction::printVoidFractionInfo() const {
  int nProcs = 0, id = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  base::MPI_Barrier();
  for (int i = 0; i < nProcs; ++i) {
    if (id == i) {
      Pout << typeName().c_str() << ":" << endl;
      for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
        Pout << "  findCellIDs[" << index << "]: " << cloud_.findCellIDs()[index] << endl;
        Pout << "  particleOverMeshNumber[" << index << "]: " << cloud_.particleOverMeshNumber()[index] << endl;
        // Pout << "  cellIDs[" << index << "]: ";
        // for (int j = 0; j < cloud_.particleOverMeshNumber()[index]; ++j) {
        //   Pout << cloud_.cellIDs()[index][j] << ", ";
        // }
        // Pout << endl;
        if (cloud_.checkFineParticle(index) || cloud_.checkMiddleParticle(index)) {
          Pout << "  voidFractions[" << index << "]: ";
          for (int j = 0; j < cloud_.particleOverMeshNumber()[index]; ++j) {
            Pout << cloud_.voidFractions()[index][j] << ", ";
          }
          Pout << endl;
          Pout << "  particleWeights[" << index << "]: ";
          for (int j = 0; j < cloud_.particleOverMeshNumber()[index]; ++j) {
            Pout << cloud_.particleWeights()[index][j] << ", ";
          }
          Pout << endl;
          Pout << "  particleVolumes[" << index << "]: ";
          for (int j = 0; j < cloud_.particleOverMeshNumber()[index]; ++j) {
            Pout << cloud_.particleVolumes()[index][j] << ", ";
          }
          Pout << endl;
        }
      }
    }
    base::MPI_Barrier();
  }  // End of procs loop
}

void mixVoidFraction::setVoidFraction() {
  // reset field
  resetVoidFraction();
  resetVolumeFraction();
  // reset particleOverMeshNumber
  base::fillTensor(cloud_.particleOverMeshNumber(), 0);
  // clear before set
  cloud_.pCloud().cellIDs().clear();
  cloud_.pCloud().voidFractions().clear();
  cloud_.pCloud().volumeFractions().clear();
  cloud_.pCloud().particleWeights().clear();
  cloud_.pCloud().particleVolumes().clear();
  // define data buffer
  // parMaps for fine and middle particles and parSets for coarse particles
  std::unordered_map<int, std::unordered_map<int, Foam::vector>> parMaps;
  std::unordered_map<int, std::unordered_set<int>> parSets;
  // set voidFraction and volumeFraction
  // 先设置大颗粒体积分数，然后设置小颗粒空隙率，因为使用高斯核函数
  // 计算小颗粒空隙率时需要使用大颗粒体积分数
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkCoarseParticle(index)) {
      // coarse particle
      parSets.insert(std::make_pair(index, std::unordered_set<int>()));
      setVolumeFractionForSingleParticle(index, parSets[index]);
    }
  }  // End loop of all coarse particles
  // 计算累计网格体积
  if (useGuassVoidFractionForMiddleParticle_) {
    getAccumulateCellVolume();
  }
  // 计算小颗粒空隙率
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    if (cloud_.checkFineParticle(index) ||
        (cloud_.checkMiddleParticle(index) && !useGuassVoidFractionForMiddleParticle_)) {
      int findCellID = cloud_.findCellIDs()[index];
      parMaps.insert(std::make_pair(index, std::unordered_map<int, Foam::vector>()));
      if (findCellID >= 0) {
        setVoidFractionForSingleParticle(index, findCellID, parMaps[index]);
      }
    } else {
      int findExpandedCellID = cloud_.findExpandedCellIDs()[index];
      parMaps.insert(std::make_pair(index, std::unordered_map<int, Foam::vector>()));
      if (findExpandedCellID >= 0) {
        setGaussVoidFractionForSingleMiddleParticle(index, findExpandedCellID, parMaps[index]);
      }
    }
  }  // End loop of all fine and middle particles
  // 在大颗粒边界处光滑空隙率场
  forAll(voidFractionNext_, cellID) {
    double volumeF = volumeFractionNext_[cellID];
    if (volumeF > Foam::SMALL && volumeF < 1.0 - Foam::SMALL) {
      voidFractionNext_[cellID] = volumeF * voidFractionNext_[cellID] + (1.0 - volumeF);
    }
  }
  voidFractionNext_.correctBoundaryConditions();
  volumeFractionNext_.correctBoundaryConditions();
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // no matter what kind of particle, need to construct default tensor
    cloud_.pCloud().cellIDs().emplace_back();
    cloud_.pCloud().voidFractions().emplace_back();
    cloud_.pCloud().volumeFractions().emplace_back();
    cloud_.pCloud().particleWeights().emplace_back();
    cloud_.pCloud().particleVolumes().emplace_back();
    if (cloud_.checkCoarseParticle(index)) {
      const auto& set = parSets[index];
      int meshNumber = set.size();
      // 检查集合中元素个数大于颗粒覆盖网格数限制
      if (meshNumber > maxCellsNumPerCoarseParticle_) {
        FatalError << __func__ << ": Big particle found " << meshNumber
                   << " cells more than permittd maximun number of cells per paticle " << maxCellsNumPerCoarseParticle_
                   << abort(FatalError);
      }
      if (meshNumber > 0) {
        // 将颗粒覆盖的当前处理器的网格数保存到 particleOverMeshNumber 中
        cloud_.particleOverMeshNumber()[index] = meshNumber;
        // allocate memory for coarse particle
        cloud_.pCloud().cellIDs()[index] = std::move(base::CITensor1(base::makeShape1(meshNumber), -1));
        cloud_.pCloud().volumeFractions()[index] = std::move(base::CDTensor1(base::makeShape1(meshNumber), 0.0));
        auto it = set.cbegin();
        int i = 0;
        for (; i < meshNumber && it != set.cend(); ++i, ++it) {
          int cellID = *it;
          // 保存颗粒覆盖的所有网格编号
          cloud_.cellIDs()[index][i] = cellID;
          // 保存 volumeFraction
          cloud_.volumeFractions()[index][i] = volumeFractionNext_[cellID];
        }
      }
    } else {
      const auto& umap = parMaps[index];
      int meshNumber = umap.size();
      if (meshNumber > 0) {
        // realloc memory according to meshNumber
        cloud_.particleOverMeshNumber()[index] = meshNumber;
        cloud_.pCloud().cellIDs()[index] = std::move(base::CITensor1(base::makeShape1(meshNumber), -1));
        cloud_.pCloud().voidFractions()[index] = std::move(base::CDTensor1(base::makeShape1(meshNumber), 0.0));
        cloud_.pCloud().particleWeights()[index] = std::move(base::CDTensor1(base::makeShape1(meshNumber), 0.0));
        cloud_.pCloud().particleVolumes()[index] = std::move(base::CDTensor1(base::makeShape1(meshNumber), 0.0));
        // 遍历 umap
        auto iter = umap.begin();
        for (int i = 0; i < meshNumber; ++i, ++iter) {
          int subCellID = iter->first;
          const Foam::vector& data = iter->second;
          // 保存颗粒覆盖的所有网格编号
          cloud_.cellIDs()[index][i] = subCellID;
          // 保存 particleWeights
          cloud_.particleWeights()[index][i] = data[0];
          // 保存 particleVolumes
          cloud_.particleVolumes()[index][i] = data[1];
          // 保存 voidFractions
          cloud_.voidFractions()[index][i] = data[1] / cloud_.mesh().V()[subCellID];
        }
      }
    }
  }  // End loop of all particles
  base::MPI_Info("mixVoidFraction: setVoidFraction - done", verbose_);
}

/*!
 * \brief 设置单个颗粒的体积分数场
 * \param index 颗粒索引
 * \param set   颗粒覆盖的网格索引的集合
 */
void mixVoidFraction::setVolumeFractionForSingleParticle(const int index, std::unordered_set<int>& set) {
  if (cloud_.checkPeriodicCells()) {
    FatalError << "Error: not support periodic check!" << abort(FatalError);
  }
  // 获取颗粒半径
  double radius = cloud_.getRadius(index);
  // 获取颗粒中心坐标
  Foam::vector particleCentre = cloud_.getPosition(index);
  // 获取到在当前 processor 上颗粒覆盖的某一个网格编号
  // 这里必须使用 findMpiCellIDs() 而不能使用 findCellIDs()
  // 因为只有颗粒中心所在的处理器的 findCellIDs()[index] >= 0
  int findMpiCellID = cloud_.findMpiCellIDs()[index];
  if (findMpiCellID >= 0) {  // particle centre is in domain
    // 获取网格中心坐标
    Foam::vector cellCentre = cloud_.mesh().C()[findMpiCellID];
    // 判断网格中心是否在颗粒中
    double fc = pointInParticle(cellCentre, particleCentre, radius);
    // 计算网格的等效半径
    double corona = 0.5 * sqrt(3.0) * cbrt(cloud_.mesh().V()[findMpiCellID]);
    // 获取网格的 corona point
    Foam::vector coronaPoint = IBVoidFraction::getCoronaPointPosition(particleCentre, cellCentre, corona);
    if (pointInParticle(coronaPoint, particleCentre, radius) < 0.0) {
      // 如果 coronaPoint 在颗粒中, 则认为整个网格在颗粒中
      volumeFractionNext_[findMpiCellID] = 0.0;
    } else {
      // 如果 coronaPoint 不在颗粒中, 则需要遍历网格的所有角点, 判断角点与网格中心是否在颗粒中
      const labelList& vertexPoints = cloud_.mesh().cellPoints()[findMpiCellID];
      double ratio = 0.125;
      double voidF = 1.0;
      // 遍历当前网格的所有角点
      forAll(vertexPoints, i) {
        if (vertexPoints[i] < 0) {
          continue;
        }
        // 获取第 i 角点坐标
        vector vertexPosition = cloud_.mesh().points()[vertexPoints[i]];
        // 判断角点是否在颗粒中
        scalar fv = pointInParticle(vertexPosition, particleCentre, radius);
        if (fc < 0.0 && fv < 0.0) {
          // 网格中心在颗粒中, 角点也在颗粒中
          voidF -= ratio;
        } else if (fc < 0.0 && fv >= 0.0) {
          // 网格中心在颗粒中, 角点不在颗粒中
          // 计算角点 vertexPosition 对体积分数的影响系数 lambda
          double lambda =
              IBVoidFraction::segmentParticleIntersection(radius, particleCentre, cellCentre, vertexPosition);
          voidF -= ratio * lambda;
        } else if (fc >= 0.0 && fv < 0.0) {
          // 网格中心不在颗粒中, 角点在颗粒中
          // 计算角点 vertexPosition 对体积分数的影响系数 lambda
          double lambda =
              IBVoidFraction::segmentParticleIntersection(radius, particleCentre, vertexPosition, cellCentre);
          voidF -= ratio * lambda;
        }
      }  // End of loop of vertexPoints
      // 保证体积分数 >= 0
      voidF = voidF < 0.0 ? 0.0 : voidF;
      voidF = voidF > 1.0 ? 1.0 : voidF;
      if (fabs(volumeFractionNext_[findMpiCellID] - 1.0) < Foam::SMALL) {
        // 如果 findMpiCellID 网格的体积分数为 1.0, 则说明第一次遍历到该网格, 可以直接赋值
        volumeFractionNext_[findMpiCellID] = voidF;
      } else {
        // 如果 findMpiCellID 网格的体积分数不为 1.0, 则说明在计算其他颗粒时候, 已经遍历到该网格
        volumeFractionNext_[findMpiCellID] -= (1.0 - voidF);
        // 保证体积分数 >= 0
        volumeFractionNext_[findMpiCellID] =
            volumeFractionNext_[findMpiCellID] < 0.0 ? 0.0 : volumeFractionNext_[findMpiCellID];
        volumeFractionNext_[findMpiCellID] =
            volumeFractionNext_[findMpiCellID] > 1.0 ? 1.0 : volumeFractionNext_[findMpiCellID];
      }
    }
    // 颗粒中心所在网格的体积分数已经计算完成, 下面开始递归构建相邻网格
    buildSetForVolumeFraction(findMpiCellID, particleCentre, radius, set);
  }  // findMpiCellID >= 0
}

/*!
 * \brief 构建颗粒覆盖的所有网格的哈希集合
 * \note 设置为递归函数,  通过哈希器将网格编号转换为哈希值, 并存入 set 中以便于搜索
 * \param cellID         <[in] 递归循环中要检索网格编号
 * \param particleCentre <[in] 颗粒中心位置
 * \param radius         <[in] 颗粒半径
 * \param set            <[in, out] 颗粒覆盖的网格索引的集合
 */
void mixVoidFraction::buildSetForVolumeFraction(const label cellID, const Foam::vector& particleCentre,
                                                const double radius, std::unordered_set<int>& set) {
  if (cellID < 0) {
    return;
  }
  set.insert(cellID);
  // 获取 cellID 网格的所有 neighbour cell 的链表
  const labelList& nc = cloud_.mesh().cellCells()[cellID];
  // 遍历 cellID 的 neighbour cell
  forAll(nc, i) {
    // 获取相邻网格索引, 以及网格中心坐标
    label neighbour = nc[i];
    if (neighbour < 0) {
      continue;
    }
    // 获取相邻网格中心坐标
    Foam::vector neighbourCentre = cloud_.mesh().C()[neighbour];
    // 判断相邻网格中心是否在颗粒中
    scalar fc = pointInParticle(neighbourCentre, particleCentre, radius);
    // 计算相邻网格的等效半径
    scalar coronaRaidus = 0.5 * sqrt(3.0) * cbrt(cloud_.mesh().V()[neighbour]);
    // 获取 corona point
    Foam::vector coronaPoint = IBVoidFraction::getCoronaPointPosition(particleCentre, neighbourCentre, coronaRaidus);
    // 如果在哈希集合中没有插入 neighbour 网格
    if (set.end() == set.find(neighbour)) {
      // 计算 neighbour 网格的体积分数
      if (pointInParticle(coronaPoint, particleCentre, radius) < 0.0) {
        // 如果相邻网格的 coronaPoint 在颗粒中, 则说明该网格完全被颗粒覆盖
        volumeFractionNext_[neighbour] = 0.0;
        // 以相邻网格为中心继续递归构建哈希集合
        buildSetForVolumeFraction(neighbour, particleCentre, radius, set);
      } else {
        // 如果相邻网格的 coronaPoint 不在颗粒中, 则需要遍历该网格的所有角点
        // 定义单个角点对空隙率的影响率
        double ratio = 0.125;
        double voidF = 1.0;
        // 获取 neighbour 网格的角点集合
        const labelList& vertexPoints = cloud_.mesh().cellPoints()[neighbour];
        /// 遍历网格 neighbour 的角点
        forAll(vertexPoints, j) {
          if (vertexPoints[j] < 0) {
            continue;
          }
          // 获取角点坐标
          Foam::vector vertexPosition = cloud_.mesh().points()[vertexPoints[j]];
          // 判断角点是否在颗粒中
          scalar fv = pointInParticle(vertexPosition, particleCentre, radius);
          if (fc < 0.0 && fv < 0.0) {  // 如果网格 neighbour 中心在颗粒中, 角点 j 也在颗粒中
            voidF -= ratio;
          } else if (fc < 0.0 && fv > 0.0) {  // 如果网格 neighbour 中心在颗粒中, 角点 j 不在颗粒中
            // 计算角点对空隙率的影响系数 lambda
            scalar lambda =
                IBVoidFraction::segmentParticleIntersection(radius, particleCentre, neighbourCentre, vertexPosition);
            voidF -= lambda * ratio;
          } else if (fc > 0.0 && fv < 0.0) {  // 如果网格 neighbour 中心不在颗粒中, 角点 j 在颗粒中
            scalar lambda =
                IBVoidFraction::segmentParticleIntersection(radius, particleCentre, vertexPosition, neighbourCentre);
            voidF -= lambda * ratio;
          }
        }  // End of loop vertexPoints
        // 保证体积分数 >= 0
        voidF = voidF < 0.0 ? 0.0 : voidF;
        voidF = voidF > 1.0 ? 1.0 : voidF;
        if (fabs(volumeFractionNext_[neighbour] - 1.0) < Foam::SMALL) {
          // 如果 neighbour 网格的体积分数为 1.0, 则说明第一次遍历到该网格, 可以直接赋值
          volumeFractionNext_[neighbour] = voidF;
        } else {
          // 如果 neighbour 网格的体积分数不为 1.0, 则说明在计算其他颗粒时候, 已经遍历到该网格
          volumeFractionNext_[neighbour] -= (1.0 - voidF);
          // 保证体积分数 >= 0
          volumeFractionNext_[neighbour] = volumeFractionNext_[neighbour] < 0.0 ? 0.0 : volumeFractionNext_[neighbour];
          volumeFractionNext_[neighbour] = volumeFractionNext_[neighbour] > 1.0 ? 1.0 : volumeFractionNext_[neighbour];
        }
        if (!(fabs(voidF - 1.0) < Foam::SMALL)) {
          // 如果体积分数不为 1.0, 则说明该 neighbour 需要递归循环构建哈希集合
          buildSetForVolumeFraction(neighbour, particleCentre, radius, set);
        }
      }
    }  // not found neighbour in hash set
  }    // End of loop neighbour cell
}

//! \brief 计算颗粒的累计网格体积
void mixVoidFraction::getAccumulateCellVolume() {
  int number = cloud_.numberOfParticles();
  // reallocate memory
  accumulateCellVolumes_ = std::move(base::CDTensor1(base::makeShape1(number), 0.0));
  Foam::vector cellPos = Foam::vector::zero;                    // 网格位置
  Foam::vector particlePos = Foam::vector::zero;                // 颗粒中心坐标
  double Vcell = 0.0;                                           // 网格体积
  double volumeF = 0.0;                                         // 大颗粒体积分数
  double GaussKernel = 0.0;                                     // 颗粒在网格处的高斯核
  base::CDTensor1 accuCellsVol(base::makeShape1(number), 0.0);  // 颗粒覆盖当前求解器的网格的累计体积
  for (int index = 0; index < cloud_.numberOfParticles(); ++index) {
    // 如果颗粒是 middle particle 同时 findExpandedCellID 位于当前求解器上
    if (cloud_.checkMiddleParticle(index) && cloud_.findExpandedCellIDs()[index] >= 0) {
      const std::unordered_set<int>& cellSet = cloud_.globalF().getExpandedCellSet(index);
      particlePos = cloud_.getPosition(index);
      // 计算累计网体积
      double accuCellVol = 0.0;
      for (const auto& cellID : cellSet) {
        // 网格 cellID 位于计算区域外部
        if (cellID < 0) {
          continue;
        }
        volumeF = volumeFractionNext_[cellID];
        // 网格 cellID 位于大颗粒内部
        if (volumeF < Foam::SMALL) {
          continue;
        }
        cellPos = cloud_.mesh().C()[cellID];
        Vcell = cloud_.mesh().V()[cellID];
        GaussKernel = GaussCore(cellPos, particlePos, GaussKernelBandWidth_);
        accuCellVol += GaussKernel * Vcell * volumeF;
      }
      accuCellsVol[index] = accuCellVol;
    }
  }
  // 主节点汇总其他节点的 accuCellsVol
  int procId = base::procId();
  int numProc = base::numProc();
  int nTag = 100;
  if (0 != procId) {
    MPI_Request rq;
    // 发送 accuCellsVol 给主节点( MPI_Isend 非阻塞)
    base::MPI_Isend(accuCellsVol, 0, nTag, &rq);
  }
  if (0 == procId) {
    std::vector<MPI_Request> rVec(numProc);
    std::vector<MPI_Status> sVec(numProc);
    std::vector<base::CDTensor1> accuCellsVolVec;
    // 接收其他节点的 accuCellsVol 信息
    for (int inode = 1; inode < numProc; ++inode) {
      accuCellsVolVec.emplace_back(base::makeShape1(number), 0);
      base::MPI_Irecv(accuCellsVolVec[inode - 1], inode, nTag, rVec.data() + inode);
    }
    // 主节点等待 Irecv 执行完成
    MPI_Waitall(numProc - 1, rVec.data() + 1, sVec.data() + 1);
    // 由主节点计算 accuCellsVol
    for (int index = 0; index < number; ++index) {
      double accuCellVol = 0.0;
      if (cloud_.checkMiddleParticle(index)) {
        for (int j = 0; j < numProc; ++j) {
          if (0 == j) {
            accuCellVol += std::max(accuCellsVol[index], 0.0);
          } else {
            accuCellVol += std::max(accuCellsVolVec[j - 1][index], 0.0);
          }
        }
        CHECK_GT(accuCellVol, 0) << __func__ << ": particle " << index << " 's accuCellVol <= 0.0";
      }
      accumulateCellVolumes_[index] = accuCellVol;
    }
  }
  // 主节点广播 accumulateCellVolume
  base::MPI_Bcast(accumulateCellVolumes_, 0);
  base::MPI_Barrier();
}

void mixVoidFraction::setGaussVoidFractionForSingleMiddleParticle(const int index, const int findExpandedCellID,
                                                                  std::unordered_map<int, Foam::vector>& parMap) {
  parMap.clear();
  double radius = cloud_.getRadius(index);                  // 颗粒半径
  Foam::vector particleCentre = cloud_.getPosition(index);  // 颗粒中心坐标
  std::unordered_set<int> expandedCellSet;                  // 扩展网格索引
  // 计算颗粒覆盖的扩展网格集合
  cloud_.voidFractionM().buildExpandedCellSet(expandedCellSet, findExpandedCellID, particleCentre, radius,
                                              GaussKernelBandWidth_ * GaussKernelScale_ / radius);
  Foam::vector cellPos = Foam::vector::zero;                    // 网格位置
  double Vcell = 0.0;                                           // 网格体积
  double volumeF = 0.0;                                         // 大颗粒体积分数
  double GaussKernel = 0.0;                                     // 颗粒在网格处的高斯核
  double effectRatio = 0.0;                                     // 颗粒影响网格系数
  double newAlpha = 1.0;                                        // 网格空隙率
  double accumulateCellVolume = accumulateCellVolumes_[index];  // 网格累计体积
  for (const int& cellID : expandedCellSet) {
    // 网格 cellID 位于计算区域外部
    if (cellID < 0) {
      continue;
    }
    volumeF = volumeFractionNext_[cellID];
    // 网格 cellID 位于大颗粒内部
    if (volumeF < Foam::SMALL) {
      continue;
    }
    cellPos = cloud_.mesh().C()[cellID];
    Vcell = cloud_.mesh().V()[cellID];
    GaussKernel = GaussCore(cellPos, particleCentre, GaussKernelBandWidth_);
    // 计算颗粒 index 对网格 cellID 的体积
    effectRatio = GaussKernel * Vcell * volumeF / accumulateCellVolume;
    // 计算空隙率
    newAlpha = voidFractionNext_[cellID] - pV(radius) * effectRatio / Vcell;
    voidFractionNext_[cellID] = newAlpha;
    // 将 cellID 插入到 map 中
    auto iter = parMap.find(cellID);
    if (parMap.end() == iter) {
      // 需要创建新的索引
      cloud_.particleOverMeshNumber()[index] += 1;
      parMap.insert(std::make_pair(cellID, Foam::vector::zero));
    }
    parMap[cellID][0] += effectRatio;               // add particleWeights
    parMap[cellID][1] += pV(radius) * effectRatio;  // add particleVolumes
    parMap[cellID][2] += pV(radius) * effectRatio;  // add particleVs
  }
}

//! \brief 设置索引为 index 的单个颗粒的空隙率
//! \note used for fine particle
void mixVoidFraction::setVoidFractionForSingleParticle(const int index, const int findCellID,
                                                       std::unordered_map<int, Foam::vector>& parMap) {
  parMap.clear();
  int numberOfMarkerPoints = dividedVoidFraction::numberOfMarkerPoints();             // 标志点个数
  int findSubCellNum = 0;                                                             // 有效标记点的个数
  double radius = cloud_.getRadius(index);                                            // 颗粒半径
  double scaleVol = weight();                                                         // 颗粒体积因子
  double scaleRadius = cbrt(porosity());                                              // 颗粒半径因子
  double avgPVol = pV(radius, scaleVol) / static_cast<double>(numberOfMarkerPoints);  // 获取平均颗粒体积（1 / 29）
  Foam::vector particleCentre = cloud_.getPosition(index);                            // 颗粒中心坐标
  Foam::vector subPosition = Foam::vector::zero;                                      // 标志点的绝对位置
  // 先计算体积, 然后再用 scaleRadius 乘以半径, 以保证在计算空隙率时候颗粒尺寸不变
  radius *= scaleRadius;
  // 遍历所有标志点，除了颗粒中心处的标志点
  for (int i = 1; i < numberOfMarkerPoints; ++i) {
    // 获取标志点的绝对位置
    subPosition = particleCentre + radius * dividedVoidFraction::offsets()[i];
    if (cloud_.checkPeriodicCells()) {
      FatalError << "Error: not support periodic check!" << abort(FatalError);
    }
    // 根据修正后的标志点坐标定位标志点所在网格的索引
    int subCellId = cloud_.locateM().findSingleCell(subPosition, findCellID);
    // 如果标志点在 domain 中, 则更新空隙率
    if (subCellId >= 0) {
      findSubCellNum += 1;
      // 获取 subCellId 网格体积
      double subCellVol = cloud_.mesh().V()[subCellId];
      double newAlpha = voidFractionNext_[subCellId] - avgPVol / subCellVol;
      // 更新空隙率场
      if (newAlpha > alphaMin_) {
        voidFractionNext_[subCellId] = newAlpha;
      } else {
        // 如果空隙率低于最低空隙率, 则直接赋值为 alphaMin_, 并累加损失的颗粒体积
        voidFractionNext_[subCellId] = alphaMin_;
        tooMuch_ += (alphaMin_ - newAlpha) * subCellVol;
      }
      auto iter = parMap.find(subCellId);
      if (parMap.end() == iter) {
        // 需要创建新的索引
        cloud_.particleOverMeshNumber()[index] += 1;
        parMap.insert(std::make_pair(subCellId, Foam::vector::zero));
      }
      parMap[subCellId][0] += 1.0 / static_cast<double>(numberOfMarkerPoints);  // add particleWeights
      parMap[subCellId][1] += avgPVol;                                          // add particleVolumes
      parMap[subCellId][2] += avgPVol;                                          // add particleVs
    }
  }  // End of marker point loop
  // 处理颗粒中心处的标志点, 将丢失的标志点对空隙率的影响作用在颗粒中心处的标志点上
  // set source for particle center
  double findCellVol = cloud_.mesh().V()[findCellID];
  double centreWeight = (numberOfMarkerPoints - findSubCellNum) * (1.0 / numberOfMarkerPoints);
  double newAlpha = voidFractionNext_[findCellID] - pV(radius, scaleVol) * centreWeight / findCellVol;
  // update voidFraction for findCellID
  if (newAlpha > alphaMin_) {
    voidFractionNext_[findCellID] = newAlpha;
  } else {
    voidFractionNext_[findCellID] = alphaMin_;
    tooMuch_ += (alphaMin_ - newAlpha) * findCellVol;
  }
  auto iter = parMap.find(findCellID);
  if (parMap.end() == iter) {
    // 需要创建新的索引
    cloud_.particleOverMeshNumber()[index] += 1;
    parMap.insert(std::make_pair(findCellID, Foam::vector::zero));
  }
  parMap[findCellID][0] += centreWeight;                         // add particleWeights
  parMap[findCellID][1] += pV(radius, scaleVol) * centreWeight;  // add particleVolumes
  parMap[findCellID][2] += pV(radius, scaleVol) * centreWeight;  // add particleVs
}

}  // namespace Foam
