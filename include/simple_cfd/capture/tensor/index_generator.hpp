#ifndef SIMPLE_CFD_CAPTURE_TENSOR_INDEX_GENERATOR_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_INDEX_GENERATOR_HPP

#include <string>
#include <simple_cfd/exception.hpp>
#include "index_variable.hpp"

namespace cfd
{

namespace detail
{

class IndexGenerator
{
private:
  char nextIndexName;

  std::string getNextIndexName()
  {
    if (nextIndexName > 'z') 
      CFD_EXCEPTION("IndexGenerator can only handle ~26 index names at present.");

    std::string result;
    result = nextIndexName++;
    return result;
  }

public:
  IndexGenerator() : nextIndexName('a')
  {
  }

  ArrayIndexVariable newArrayIndexVariable(const std::size_t limit)
  {
    return ArrayIndexVariable(getNextIndexName(), limit);
  }

  TensorIndexVariable newTensorIndexVariable(const std::size_t limit)
  {
    return TensorIndexVariable(getNextIndexName(), limit);
  }

};

}

}

#endif
