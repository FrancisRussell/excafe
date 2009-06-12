#ifndef SIMPLE_CFD_FE_BINARY_FUNCTION_HPP
#define SIMPLE_CFD_FE_BINARY_FUNCTION_HPP

#include <utility>
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
  virtual double evaluate(const mesh<TriangularCell>& m, const MeshEntity& entity, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const = 0;
};

template<typename T>
class ScaledFEBinaryFunction : public FEBinaryFunction<typename T::cell_type>
{
public:
  typedef T ScaledFunctionType;
  typedef typename ScaledFunctionType::cell_type cell_type;
  typedef typename ScaledFunctionType::vertex_type vertex_type;
  typedef typename ScaledFunctionType::finite_element_t finite_element_t;

  const ScaledFunctionType function;
  const double scaleFactor;

  ScaledFEBinaryFunction(const ScaledFunctionType& f, const double _scaleFactor) : function(f), scaleFactor(_scaleFactor)
  {
  }

  virtual const finite_element_t* getTestFunction() const
  {
    return function.getTestFunction();
  }

  virtual const finite_element_t* getTrialFunction() const
  {
    return function.getTrialFunction();
  }

  virtual double evaluate(const mesh<TriangularCell>& m, const MeshEntity& entity, const std::size_t testDof, const std::size_t trialDof, const vertex_type& location) const
  {
    return function.evaluate(m, entity, testDof, trialDof, location) * scaleFactor;
  }
};

}
#endif
