#include "./cfdem_tools.h"
#include "sub_model/force_model/force_model.h"

namespace Foam {

/*!
 * \brief 根据 Zhou et al. 2010,JFM 的论文, 定义了三种模型, setI, setII, setIII
 *   "Bfull" 模型对应 setI
 *   "A" 模型对应 setII
 *   "B" 模型对应 setIII
 */
void cfdemTools::checkModelType(const cfdemCloud& cloud) {
  const std::string& modelType = cloud.modelType();
  if ("Bfull" == modelType || "A" == modelType) {
    // // 必须使用压力梯度力
    // if (!isUsedForceModel(cloud, "gradPForce") && !isUsedForceModel(cloud, "mixGradPForce")) {
    //   FatalError << __func__ << ": gradPForce or mixGradPForce not found with model type " << modelType
    //              << abort(FatalError);
    // }
    // // 必须使用粘性力
    // if (!isUsedForceModel(cloud, "viscForce") && !isUsedForceModel(cloud, "mixViscForce")) {
    //   FatalError << __func__ << ": viscForce or mixViscForce not found with model type " << modelType
    //              << abort(FatalError);
    // }
  } else if (modelType == "B") {
    // 必须使用 Archimedes force, 且只能指定 Archimedes or mixArchimedes 中任意一个
    if (isUsedForceModel(cloud, "Archimedes") && isUsedForceModel(cloud, "mixArchimedes")) {
      FatalError << __func__ << ": Archimedes or mixArchimedes model both found! You must use only one of them\n"
                 << abort(FatalError);
    } else if ((!isUsedForceModel(cloud, "Archimedes")) && (!isUsedForceModel(cloud, "mixArchimedes"))) {
      FatalError << __func__ << ": Archimedes and mixArchimedes model both not found! You must use only one of them\n"
                 << abort(FatalError);
    }
    // 不能使用压力梯度力
    if (isUsedForceModel(cloud, "gradPForce")) {
      FatalError << __func__ << ": can not use gradPForce with model type B!\n" << abort(FatalError);
    }
    // 不能使用粘性力
    if (isUsedForceModel(cloud, "viscForce")) {
      FatalError << __func__ << ": can not use viscForce with model type B!\n" << abort(FatalError);
    }
  } else if (modelType == "none") {
    Info << "\nsolving volume averaged Navier Stokes equations of type none\n" << endl;
  } else {
    FatalError << "\nno suitable model type specified:" << modelType << endl << abort(FatalError);
  }
}

bool cfdemTools::isUsedForceModel(const cfdemCloud& cloud, const std::string& forceModelName) {
  bool isUsed = false;
  for (const auto& sPtr : cloud.forceModels()) {
    if (forceModelName == sPtr->typeName()) {
      isUsed = true;
      break;
    }
  }
  return isUsed;
}

}  // namespace Foam
