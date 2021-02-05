#ifndef __AVERAGING_MODEL_INL_H__
#define __AVERAGING_MODEL_INL_H__

#include "./averaging_model.h"

namespace Foam {

/*!
 * \brief 计算局部累加矢量场
 * \param valueField   <[in, out] 需要被局部累加的场
 * \param value        <[in] 用于局部累加的颗粒数据(lagrange value)
 * \param weight       <[in] 用于局部累加的权重系数(lagrange value)
 */
template <int nDim, typename DType, typename Device, typename Alloc>
void averagingModel::setVectorFieldSum(volVectorField& valueField,
                                       const base::Tensor<nDim, DType, Device, Alloc>& value,
                                       const std::vector<base::CDTensor1>& weight) {
  CHECK(2 == nDim) << __func__ << ": error tensor dimension";
  CHECK_EQ(value.size(1), 3) << __func__ << ": vector field's dimension must equal to 3";
  for (int index = 0; index < cloud_.numberOfParticles(); index++) {
    if (cloud_.checkFineParticle(index) || cloud_.checkMiddleParticle(index)) {
      for (int subCell = 0; subCell < cloud_.particleOverMeshNumber()[index]; subCell++) {
        int cellID = cloud_.cellIDs()[index][subCell];
        if (cellID >= 0) {
          // get value vector
          Foam::vector valueVec(value[index][0], value[index][1], value[index][2]);
          valueField[cellID] += valueVec * weight[index][subCell];
        }
      }
    }
  }
}

}  // namespace Foam

#endif  // __AVERAGING_MODEL_INL_H__
