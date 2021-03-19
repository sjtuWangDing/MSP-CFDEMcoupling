#ifndef __CLOUD_POST_PROC__
#define __CLOUD_POST_PROC__

#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include "base/logging.h"
#include "base/string_utils.h"
#include "vector.H"

struct Atom {
  Atom() : id_(-1), radius_(0.0), pos_(Foam::vector::zero), vel_(Foam::vector::zero) {}
  Atom(int id, double radius, Foam::vector pos, Foam::vector vel) : id_(id), radius_(radius), pos_(pos), vel_(vel) {}
  Atom(const std::string& dataStr) {
    std::vector<std::string> subStrVec = base::split(dataStr, " ", 18);
    id_ = std::stoi(subStrVec[0]);
    radius_ = std::stod(subStrVec[17]);
    pos_ = Foam::vector(std::stod(subStrVec[2]), std::stod(subStrVec[3]), std::stod(subStrVec[4]));
    vel_ = Foam::vector(std::stod(subStrVec[8]), std::stod(subStrVec[9]), std::stod(subStrVec[10]));
  }
  friend std::ostream& operator<<(std::ostream& os, const Atom& atom) {
    os << atom.id_ << " " << atom.radius_ << " ";
    os << "[" << atom.pos_[0] << ", " << atom.pos_[1] << ", " << atom.pos_[2] << "] ";
    os << "[" << atom.vel_[0] << ", " << atom.vel_[1] << ", " << atom.vel_[2] << "]";
    return os;
  }
  int id_;
  double radius_;
  Foam::vector pos_;
  Foam::vector vel_;
};

struct Node {
  Node() : timeStep_(0.0), atomNum_(0) {}
  double timeStep_;
  int atomNum_;
  std::vector<Atom> atomVec_;
};

class cfdemPostProc {
 public:
  //! \brief 从指定配置文件中加载配置信息
  std::vector<Node> ReadDmpRunFile(const std::string& dmpFileName, const double DEMDt) {
    // 打开文件
    std::ifstream ifs(dmpFileName.c_str());
    CHECK(ifs) << "load dump file " << dmpFileName.c_str() << " failed!" << std::endl;

    // dump file string
    std::string dmpLine;
    bool isTimeStep = false, isNumber = false, isAtoms = false;
    std::vector<Node> NodeVec;
    double startTimeStep = 0;
    // 默认以'\n'为间隔依次读取文件流的一行
    for (int lineNo = 0; getline(ifs, dmpLine, '\n'); ++lineNo) {
      dmpLine = base::trim(dmpLine);
      std::vector<std::string> subStrVec = base::split(dmpLine, ":", 1);
      if (subStrVec.size() == 2) {
        isTimeStep = false;
        isNumber = false;
        isAtoms = false;
        if (!strcasecmp(subStrVec[1].c_str(), "TIMESTEP")) {
          isTimeStep = true;
          NodeVec.emplace_back();
        } else if (!strcasecmp(subStrVec[1].c_str(), "NUMBER OF ATOMS")) {
          isNumber = true;
        } else if (!strcasecmp(subStrVec[1].c_str(),
                               "ATOMS id type x y z ix iy iz vx vy vz fx fy fz omegax omegay omegaz radius")) {
          isAtoms = true;
        } else {
          isTimeStep = false;
          isNumber = false;
          isAtoms = false;
        }
        continue;
      }
      auto iter = --NodeVec.end();
      if (isTimeStep) {
        if (1 == NodeVec.size()) {
          startTimeStep = std::stod(dmpLine);
        }
        iter->timeStep_ = (std::stod(dmpLine) - startTimeStep) * DEMDt;
        std::cout << startTimeStep << ", " << iter->timeStep_ << std::endl;
      }
      if (isNumber) {
        iter->atomNum_ = std::stoi(dmpLine);
      }
      if (isAtoms) {
        iter->atomVec_.emplace_back(dmpLine);
      }
    }
    // 关闭dmp文件
    ifs.close();
    return NodeVec;
  }

  std::vector<std::pair<double, Foam::vector>> volAverage(Foam::vector recStartPoint, Foam::vector recDims,
                                                          double startTime, const std::vector<Node>& NodeVec) {
    std::vector<std::pair<double, Foam::vector>> volVec;
    // 遍历时间步
    for (const auto& node : NodeVec) {
      if (node.timeStep_ < startTime) {
        continue;
      }
      Foam::vector sumVolVel = Foam::vector::zero;
      double sumVol = 0;
      for (const Atom& atom : node.atomVec_) {
        if (checkParticleInRec(atom.pos_, recStartPoint, recDims)) {
          sumVolVel += pV(atom.radius_) * atom.vel_;
          sumVol += pV(atom.radius_);
        }
      }
      if (sumVol > Foam::SMALL) {
        volVec.emplace_back(std::make_pair(node.timeStep_, sumVolVel / sumVol));
        // std::cout << volVec[volVec.size() - 1].first << " " << volVec[volVec.size() - 1].second[0] << " "
        //           << volVec[volVec.size() - 1].second[1] << " " << volVec[volVec.size() - 1].second[2] << std::endl;
      }
    }
    return volVec;
  }

  Foam::vector timeAverage(Foam::vector recStartPoint, Foam::vector recDims, const double startTime, const double DEMDt,
                           const std::vector<Node>& NodeVec) {
    std::vector<std::pair<double, Foam::vector>> volVec;
    Foam::vector timeAveVel = Foam::vector::zero;
    Foam::vector timeSumVel = Foam::vector::zero;
    double timeSum = 0.0;
    // 遍历时间步
    for (const auto& node : NodeVec) {
      if (node.timeStep_ < startTime) {
        continue;
      }
      // 计算 Rect 中体积平均速度
      Foam::vector sumVolVel = Foam::vector::zero;
      double sumVol = 0;
      for (const Atom& atom : node.atomVec_) {
        if (checkParticleInRec(atom.pos_, recStartPoint, recDims)) {
          Foam::vector atomVel(atom.vel_[0], atom.vel_[1], atom.vel_[2] > 0.21 ? 0.21 : atom.vel_[2]);
          sumVolVel += pV(atom.radius_) * atomVel;
          sumVol += pV(atom.radius_);
        }
      }
      if (sumVol > Foam::SMALL) {
        // 累加时间平均速度
        timeSumVel += (sumVolVel / sumVol) * 100 * DEMDt;
        // 累加时间
        timeSum += 100 * DEMDt;
      }
    }
    if (timeSum > Foam::SMALL) {
      timeAveVel = timeSumVel / timeSum;
    }
    return timeAveVel;
  }

  Foam::vector timeAverage(Foam::vector cycStartPoint, const double radius, const double dimZ, const double startTime,
                           const double DEMDt, const std::vector<Node>& NodeVec) {
    std::vector<std::pair<double, Foam::vector>> volVec;
    Foam::vector timeAveVel = Foam::vector::zero;
    Foam::vector timeSumVel = Foam::vector::zero;
    double timeSum = 0.0;
    // 遍历时间步
    for (const auto& node : NodeVec) {
      if (node.timeStep_ < startTime) {
        continue;
      }
      // 计算 Rect 中体积平均速度
      Foam::vector sumVolVel = Foam::vector::zero;
      double sumVol = 0;
      for (const Atom& atom : node.atomVec_) {
        if (checkParticleInCyc(atom.pos_, cycStartPoint, radius, dimZ)) {
          Foam::vector atomVel(atom.vel_[0], atom.vel_[1], atom.vel_[2] > 0.21 ? 0.21 : atom.vel_[2]);
          sumVolVel += pV(atom.radius_) * atomVel;
          sumVol += pV(atom.radius_);
        }
      }
      if (sumVol > Foam::SMALL) {
        // 累加时间平均速度
        timeSumVel += (sumVolVel / sumVol) * 100 * DEMDt;
        // 累加时间
        timeSum += 100 * DEMDt;
      }
    }
    if (timeSum > Foam::SMALL) {
      timeAveVel = timeSumVel / timeSum;
    }
    return timeAveVel;
  }

  inline double pV(const double& radius) { return 4.0 * M_PI * radius * radius * radius / 3.0; }

  inline bool checkParticleInRec(const Foam::vector& particlePos, const Foam::vector& recStartPoint,
                                 const Foam::vector& recDims) {
    return (particlePos[0] >= recStartPoint[0] && particlePos[1] >= recStartPoint[1] &&
            particlePos[2] >= recStartPoint[2]) &&
           ((particlePos[0] - recStartPoint[0]) <= recDims[0] && (particlePos[1] - recStartPoint[1]) <= recDims[1] &&
            (particlePos[2] - recStartPoint[2]) <= recDims[2]);
    //  (fabs(particlePos[0] - recStartPoint[0]) <= recDims[0] && fabs(particlePos[1] - recStartPoint[1]) <= recDims[1]
    //  &&
    //   fabs(particlePos[2] - recStartPoint[2]) <= recDims[2]);
  }

  inline bool checkParticleInCyc(const Foam::vector& particlePos, const Foam::vector& cycStartPoint,
                                 const double radius, const double dimZ) {
    double posX = particlePos.x() - cycStartPoint.x();
    double posY = particlePos.y() - cycStartPoint.y();
    double posZ = particlePos.z() - cycStartPoint.z();
    return posZ >= 0 && posZ <= dimZ && pow(posX * posX + posY * posY, 0.5) < radius + Foam::SMALL;
  }
};

#endif  // __CLOUD_POST_PROC__
