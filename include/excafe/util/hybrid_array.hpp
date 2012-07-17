#ifndef EXCAFE_UTIL_HYBRID_ARRAY_HPP
#define EXCAFE_UTIL_HYBRID_ARRAY_HPP

#include <cstddef>

namespace excafe
{

namespace util
{

template<typename T, std::size_t N>
class HybridArray
{
public:
  typedef T value_type;
  static const std::size_t max_stack_size = N;

private:
   HybridArray(const HybridArray&);
   HybridArray& operator=(const HybridArray&);

   value_type stackData[max_stack_size];
   value_type* const dataPointer;

public:
   HybridArray(const std::size_t size) : 
     dataPointer(size > max_stack_size ? new value_type[size] : stackData)
   {
   }

   const value_type& operator[](const std::size_t i) const
   {
     return dataPointer[i];
   }

   value_type& operator[](const std::size_t i)
   {
     return dataPointer[i];
   }

   value_type* get()
   {
     return dataPointer;
   }

   const value_type* get() const
   {
     return dataPointer;
   }

   ~HybridArray()
   {
     if (dataPointer != stackData)
       delete[] dataPointer;
   }
};

}

}

#endif
