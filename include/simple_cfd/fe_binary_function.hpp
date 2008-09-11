#ifndef SIMPLE_CFD_FE_BINARY_FUNCTION_HPP
#define SIMPLE_CFD_FE_BINARY_FUNCTION_HPP

#include "finite_element.hpp"

namespace cfd
{

template<typename C>
class FEBinaryFunction
{
private:
  typedef C cell_type;
  typedef typename cell_type::vertex_type vertex_type;
  typedef finite_element<cell_type> finite_element_t;

public:
  virtual finite_element_t* getTestFunction() const;
  virtual finite_element_t* getTrialFunction() const;
  virtual double evaluate(const cell_type& cell, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const;
};

}
#endif
