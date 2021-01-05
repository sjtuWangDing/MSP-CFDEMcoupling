#ifndef __GLOBAL_FORCE_H__
#define __GLOBAL_FORCE_H__

#include "cloud/cfdem_cloud.h"

namespace Foam {

class globalForce {
 public:
  globalForce(cfdemCloud& cloud);

  ~globalForce();

  inline void resetImpParticleForce() {
    impParticleForce_ == dimensionedVector("zero", impParticleForce_.dimensions(), vector::zero);
  }

  inline void resetExpParticleForcee() {
    expParticleForce_ == dimensionedVector("zero", expParticleForce_.dimensions(), vector::zero);
  }

  inline const volVectorField& impParticleForce() const { return impParticleForce_; }

  inline const volVectorField& expParticleForce() const { return expParticleForce_; }

  inline volVectorField& impParticleForce() { return impParticleForce_; }

  inline volVectorField& expParticleForce() { return expParticleForce_; }

 private:
  cfdemCloud& cloud_;

  //! \brief 颗粒隐式力的总和 [N]
  volVectorField impParticleForce_;

  //! \brief 颗粒显式力的总和 [N]
  volVectorField expParticleForce_;
};

}  // namespace Foam

#endif  // __GLOBAL_FORCE_H__
