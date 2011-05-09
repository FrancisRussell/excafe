#ifndef SIMPLE_CFD_CODEGEN_UFL_EXPRESSION_PROVIDER_HPP
#define SIMPLE_CFD_CODEGEN_UFL_EXPRESSION_PROVIDER_HPP

#include <string>
#include <simple_cfd/cse/expression_provider.hpp>
#include <simple_cfd/capture/assembly/scalar_placeholder.hpp>

namespace cfd
{

namespace codegen
{

class UFLExpressionProvider : public cse::ExpressionProvider<detail::ScalarPlaceholder>
{
  virtual std::string getLValue(const unsigned index)
  {
    return "lvalue";
  }
  
  std::string getRValue(const variable_t& var)
  {
    return "rvalue";
  }
};

}

}

#endif
