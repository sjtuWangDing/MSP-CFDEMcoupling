#include "./global_force.h"

namespace Foam {

globalForce::globalForce(cfdemCloud& cloud)
    : cloud_(cloud),
      impParticleForce_(
          IOobject("impParticleForce", cloud.mesh().time().timeName(), cloud.mesh(), IOobject::READ_IF_PRESENT,
                   IOobject::AUTO_WRITE),
          cloud.mesh(),
          dimensionedVector("zero", dimensionSet(1, 1, -2, 0, 0), vector(0, 0, 0))  // [N] == [kg * m / s^2]
          ),
      expParticleForce_(
          IOobject("expParticleForce", cloud.mesh().time().timeName(), cloud.mesh(), IOobject::READ_IF_PRESENT,
                   IOobject::AUTO_WRITE),
          cloud.mesh(),
          dimensionedVector("zero", dimensionSet(1, 1, -2, 0, 0), vector(0, 0, 0))  // [N] == [kg * m / s^2]
      ) {}

globalForce::~globalForce() {}

}  // namespace Foam
