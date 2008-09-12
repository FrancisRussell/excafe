#ifndef SIMPLE_CFD_FE_BINARY_FUNCTION_HPP
#define SIMPLE_CFD_FE_BINARY_FUNCTION_HPP

#include "finite_element.hpp"

namespace cfd
{

template<typename C>
class FEBinaryFunction
{
public:
  typedef C cell_type;
  typedef typename cell_type::vertex_type vertex_type;
  typedef finite_element<cell_type> finite_element_t;

  virtual const finite_element_t* getTestFunction() const = 0;
  virtual const finite_element_t* getTrialFunction() const = 0;
  virtual double evaluate(const cell_type& cell, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const = 0;
};

}
#endif
