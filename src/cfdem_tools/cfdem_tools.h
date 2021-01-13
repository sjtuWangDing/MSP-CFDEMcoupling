#ifndef __CFDEM_TOOLS_H__
#define __CFDEM_TOOLS_H__

#include "cloud/cfdem_cloud.h"

namespace Foam {

class cfdemTools {
 public:
  static void checkModelType(const cfdemCloud& cloud);

 protected:
  static bool isUsedForceModel(const cfdemCloud& cloud, const std::string& forceModelName);
};

}  // namespace Foam

#endif  // __CFDEM_TOOLS_H__
