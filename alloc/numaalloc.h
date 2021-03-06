
#include <numa.h>

template <class T>
class NumaAlloc {
  const int node;

 public:
  typedef T value_type;
  typedef std::true_type propagate_on_container_move_assignment;
  typedef std::true_type is_always_equal;

  NumaAlloc(int node = 0) noexcept : node{node} {};

  T* allocate(size_t num) {
    auto ret = numa_alloc_onnode(num * sizeof(T), node);
    if (!ret) throw std::bad_alloc();
    return reinterpret_cast<T*>(ret);
  }

  void deallocate(T* p, size_t num) noexcept { numa_free(p, num * sizeof(T)); }
};

template <class T1, class T2>
bool operator==(const NumaAlloc<T1>&, const NumaAlloc<T2>&) noexcept {
  return true;
}

template <class T1, class T2>
bool operator!=(const NumaAlloc<T1>&, const NumaAlloc<T2>&) noexcept {
  return false;
}
