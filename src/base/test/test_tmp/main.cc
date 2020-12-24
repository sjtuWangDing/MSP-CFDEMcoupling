#include "base/memory/tmp.h"
using std::cout;
using std::endl; 

class BG : public base::RefCounter {
public:
  double data_[1000000];
};

// 编译器默认会RVO
base::Tmp<BG> test01() {
  base::Tmp<BG> tmp1(new BG);
  return tmp1;
}

// 强制编译器调用移动构造函数
base::Tmp<BG> test02() {
  base::Tmp<BG> tmp1(new BG);
  return std::move(tmp1);
}

int main() {
  base::Tmp<BG> tmp1 = test01();
  base::Tmp<BG> tmp2 = test02();
  return 0;
}
