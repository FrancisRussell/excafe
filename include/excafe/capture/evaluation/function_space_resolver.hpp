#ifndef EXCAFE_CAPTURE_EVALUATION_FUNCTION_SPACE_RESOLVER_HPP
#define EXCAFE_CAPTURE_EVALUATION_FUNCTION_SPACE_RESOLVER_HPP

#include <map>
#include <cstddef>
#include <utility>
#include <excafe/capture/capture_fwd.hpp>
#include <excafe/capture/fields/function_space_expr.hpp>
#include <excafe/capture/fields/function_space_addition.hpp>
#include <excafe/capture/fields/function_space_visitor.hpp>
#include <excafe/capture/fields/function_space_mesh.hpp>
#include <excafe/dof_map_builder.hpp>
#include <excafe/exception.hpp>

namespace excafe
{

namespace detail
{

template<std::size_t D>
class FunctionSpaceResolver : public FunctionSpaceVisitor
{
private:
  static const std::size_t dimension = D;
  typedef detail::FunctionSpaceExpr* function_space_ptr;

  Scenario<dimension>& scenario;
  std::map< function_space_ptr, DofMap<dimension> >& functionSpaceMap;

  bool alreadyResolved(detail::FunctionSpaceExpr& f) const
  {
    return functionSpaceMap.find(&f) != functionSpaceMap.end();
  }

  DofMap<dimension> getDofMap(FunctionSpaceExpr& e) const
  {
    typename std::map< function_space_ptr, DofMap<dimension> >::const_iterator mapIter = functionSpaceMap.find(&e);
    assert(mapIter != functionSpaceMap.end());
    return mapIter->second;
  }

  void setDofMap(FunctionSpaceExpr& e, const DofMap<dimension>& m)
  {
    assert(functionSpaceMap.find(&e) == functionSpaceMap.end());
    functionSpaceMap.insert(std::make_pair(&e, m));
  }

public:
  FunctionSpaceResolver(Scenario<dimension>& _scenario, std::map< function_space_ptr, DofMap<dimension> >& _functionSpaceMap) :
    scenario(_scenario), functionSpaceMap(_functionSpaceMap)
  {
  }

  virtual void visit(FunctionSpaceAddition& f)
  {
    if (alreadyResolved(f))
      return;

    f.getLeft().accept(*this);
    f.getRight().accept(*this);

    setDofMap(f, getDofMap(f.getLeft()) + getDofMap(f.getRight()));
  }

  virtual void visit(FunctionSpaceMesh& f)
  {
    if (alreadyResolved(f))
      return;

    DofMapBuilder<dimension> mapBuilder(scenario.getMesh());
    mapBuilder.addFiniteElement(scenario.getElement(f.getElement()));
    setDofMap(f, mapBuilder.getDofMap());
  }

  virtual void visit(FunctionSpaceUndefined& f)
  {
    CFD_EXCEPTION("Tried to evaluate undefined function space");
  }
};

}

}

#endif
