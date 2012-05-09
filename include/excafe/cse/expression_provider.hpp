#ifndef EXCAFE_CSE_EXPRESSION_PROVIDER_HPP
#define EXCAFE_CSE_EXPRESSION_PROVIDER_HPP

#include <string>

namespace excafe
{

namespace cse
{

template<typename V>
class ExpressionProvider
{
public:
  typedef V variable_t;

  virtual std::string getLValue(const unsigned index) = 0;
  virtual std::string getRValue(const variable_t& var) = 0;
};

}

}

#endif
