#include <iostream>
#include <memory>
#include <vector>
#include "base/tensor/tensor.h"
using std::cout;
using std::endl;

// test 1-dim tensor
void test01() {
  base::CITensor1 t1(base::makeShape1(20));
  cout << (void*)t1.dptr_ << endl;
  cout << (void*)t1.optr_ << endl;
  cout << t1.mSize() << endl;
  cout << t1.stride_ << endl;
  cout << t1.memSize<0>() << endl;
  for (int i = 0; i < t1.mSize(); ++i) {
    t1.dptr_[i] = 10;
  }
}

// test 2-dim tensor
void test02() {
  base::CITensor2 t(base::makeShape2(10, 20));
  cout << (void*)t.dptr_ << endl;
  cout << (void*)t.optr_ << endl;
  cout << t.mSize() << endl;
  cout << t.stride_ << endl;
  cout << t.memSize<0>() << endl;
  cout << t.memSize<1>() << endl;
  std::fill_n(t.dptr_, t.mSize(), 5);
  for (int i = 0; i < t.mSize(); ++i) {
    cout << t.dptr_[i] << ", ";
  }
  cout << endl;
  auto t1 = t[0];
  cout << (void*)t1.dptr_ << endl;
  cout << (void*)t1.optr_ << endl;
  cout << t1.mSize() << endl;
  for (int i = 0; i < t.size(0); ++i) {
    for (int j = 0; j < t.size(1); ++j) {
      cout << t[i][j] << ", ";
    }
    cout << endl;
  }
  cout << endl;
}

// test move construct
void test03() {
  base::CITensor1 t2;
  base::CITensor1 t1(base::makeShape1(20));
  t2 = std::move(t1);
  cout << (void*)t1.dptr_ << endl;
  cout << (void*)t1.optr_ << endl;
  cout << (void*)t2.dptr_ << endl;
  cout << (void*)t2.optr_ << endl;
}

// test stl container
void test04() {
  // std::vector<std::unique_ptr<double>> v1(10, std::unique_ptr<double>(new double));
  std::vector<std::unique_ptr<double>> v1;
  v1.push_back(std::unique_ptr<double>(new double(10)));
  v1.push_back(std::unique_ptr<double>(new double(20)));
  v1.push_back(std::unique_ptr<double>(new double(30)));
  for (auto& ptr : v1) {
    cout << *ptr << endl;
  }

  std::vector<base::CITensor1> v2;
  base::CITensor1 t1(base::makeShape1(2));
  base::CITensor1 t2(base::makeShape1(3));
  base::CITensor1 t3(base::makeShape1(4));
  std::fill_n(t1.dptr_, t1.mSize(), 10);
  std::fill_n(t2.dptr_, t2.mSize(), 20);
  std::fill_n(t3.dptr_, t3.mSize(), 30);
  cout << (void*)t1.dptr_ << endl;
  cout << (void*)t2.dptr_ << endl;
  cout << (void*)t3.dptr_ << endl;
  v2.push_back(std::move(t1));
  v2.push_back(std::move(t2));
  v2.push_back(std::move(t3));
  for (const auto& tensor : v2) {
    for (int i = 0; i < tensor.mSize(); ++i) {
      cout << tensor[i] << ", ";
    }
    cout << endl;
  }
  cout << endl;
}

void reallocateTensorVector() {
  std::vector<base::CITensor1> v1;
  int i = 0;
  int particleNumber = 0;
  int meshNumber = 100;
  while (i++ < 1000000) {
    particleNumber += 2;
    meshNumber += 100;
    v1.clear();
    for (int index = 0; index < particleNumber; ++index) {
      v1.emplace_back(base::makeShape1(meshNumber), -1);
    }
    cout << "memory size: " << (meshNumber * particleNumber * sizeof(int)) / (1024 * 1024) << "MB" << endl;
  }
}

int main(int argc, const char* argv[]) {
  // test01();
  // test02();
  // test03();
  // test04();
  reallocateTensorVector();
  return 0;
}
