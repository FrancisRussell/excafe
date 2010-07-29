#ifndef SIMPLE_CFD_UTIL_HYBRID_ARRAY_HPP
#define SIMPLE_CFD_UTIL_HYBRID_ARRAY_HPP

#include <cstddef>
#include <boost/utility.hpp>
#include <boost/scoped_array.hpp>

namespace cfd
{

namespace util
{

template<typename T, std::size_t N>
class HybridArray : boost::noncopyable
{
public:
  typedef T value_type;
  static const std::size_t max_stack_size = N;

private:
   value_type stackAllocated[max_stack_size];
   const boost::scoped_array<value_type> heapAllocated;
   value_type* const dataPointer;

   static value_type* allocateDynamic(const std::size_t size)
   {
     if (size > max_stack_size)
       return new value_type[size];
     else
       return NULL;
   }

   value_type* getData()
   {
     return heapAllocated.get() != NULL ? heapAllocated.get() : stackAllocated;
   }

public:
   HybridArray(const std::size_t size) : heapAllocated(allocateDynamic(size)), dataPointer(getData())
   {
   }

   const value_type& operator[](const std::ptrdiff_t i) const
   {
     return dataPointer[i];
   }

   value_type& operator[](const std::ptrdiff_t i)
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
};

}

}

#endif
