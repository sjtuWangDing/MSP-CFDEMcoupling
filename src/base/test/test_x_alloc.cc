#include <iostream>
#include "../memory/x_alloc.h"
using namespace std;

union Obj {
  union Obj* next;
  char data[1];
};

void print(void* memory, size_t size) {
  for (size_t i = 0; i < size; ++i) {
    cout << ((char*)memory)[i] << ", ";
  }
  cout << endl;
}

void test01() {
  void* memory1 = malloc(16);
  void* memory2 = malloc(16);
  std::fill_n((char*)memory1, 16, '0');
  std::fill_n((char*)memory2, 16, '0');
  ((Obj*)memory1)->next = (Obj*)memory2;
  ((Obj*)memory2)->next = (Obj*)0;
  cout << *((unsigned long*)memory1) << endl;
  cout << (unsigned long)memory2 << endl;
  cout << (void*)memory1 << endl;
  cout << (void*)((Obj*)memory1)->data << endl;
  free(memory1);
  free(memory2);
  memory1 = 0;
  memory2 = 0;
}

void test02() {
  base::simple_alloc<int, base::malloc_alloc> alloc1;
  int* p = alloc1.allocate(20);
  alloc1.deallocate(p);
  p = 0;
}

void test03() {
  base::simple_alloc<int, base::single_thread_alloc> alloc1;
  int* p = alloc1.allocate(12);
  alloc1.deallocate(p, 12);
  p = 0;
}

int main(int argc, const char * argv[]) {
  test01();
  test02();
  test03();
  return 0;
}
