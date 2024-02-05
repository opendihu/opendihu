#pragma once

#include <cstdlib>

// CC-BY-SA by user1071136, source: https://stackoverflow.com/a/8545389/10290071
template <typename T, std::size_t N = 16> class AlignmentAllocator {
public:
  typedef T value_type;
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;

  typedef T *pointer;
  typedef const T *const_pointer;

  typedef T &reference;
  typedef const T &const_reference;

public:
  inline AlignmentAllocator() throw() {}

  template <typename T2>
  inline AlignmentAllocator(const AlignmentAllocator<T2, N> &) throw() {}

  inline ~AlignmentAllocator() throw() {}

  inline pointer adress(reference r) { return &r; }

  inline const_pointer adress(const_reference r) const { return &r; }

  inline pointer allocate(size_type n) {
    return (pointer)aligned_alloc(N, n * sizeof(value_type));
  }

  inline void deallocate(pointer p, size_type) { free(p); }

  inline void construct(pointer p, const value_type &wert) {
    new (p) value_type(wert);
  }

  inline void destroy(pointer p) { p->~value_type(); }

  inline size_type max_size() const throw() {
    return size_type(-1) / sizeof(value_type);
  }

  template <typename T2> struct rebind {
    typedef AlignmentAllocator<T2, N> other;
  };

  bool operator!=(const AlignmentAllocator<T, N> &other) const {
    return !(*this == other);
  }

  // Returns true if and only if storage allocated from *this
  // can be deallocated from other, and vice versa.
  // Always returns true for stateless allocators.
  bool operator==(const AlignmentAllocator<T, N> &other) const { return true; }
};
