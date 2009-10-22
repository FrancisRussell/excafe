#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_FUNCTION_SPACE_RESOLVER_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_FUNCTION_SPACE_RESOLVER_HPP

#include <map>
#include <cstddef>
#include <utility>
#include <simple_cfd/capture/fields/function_space_expr.hpp>
#include <simple_cfd/capture/fields/function_space_addition.hpp>
#include <simple_cfd/capture/fields/function_space_visitor.hpp>
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace detail
{

template<std::size_t D>
class FunctionSpaceResolver : public FunctionSpaceVisitor
{
private:
  static const std::size_t dimension = D;
  typedef detail::FunctionSpaceExpr* function_space_ptr;

  std::map< function_space_ptr, DofMap<dimension> >& functionSpaceMap;

  bool alreadyResolved(const detail::FunctionSpaceExpr& f) const
  {
    return functionSpaceMap.find(&f) != functionSpaceMap.end();
  }

  DofMap<dimension> getDofMap(const FunctionSpaceExpr& e) const
  {
    typename std::map< function_space_ptr, DofMap<dimension> >::const_iterator mapIter = functionSpaceMap.find(&e);
    assert(mapIter != functionSpaceMap.end());
    return mapIter->second;
  }

  void setDofMap(FunctionSpaceExpr& e, const DofMap<dimension>& m)
  {
    functionSpaceMap.insert(std::make_pair(&e, m));
  }

public:
  FunctionSpaceResolver(std::map< function_space_ptr, DofMap<dimension> >& _functionSpaceMap) :
    functionSpaceMap(_functionSpaceMap)
  {
  }

  virtual void visit(FunctionSpaceAddition& f)
  {
    if (alreadyResolved(f))
      return;

    setDofMap(f, getDofMap(f.getLeft()) + getDofMap(f.getRight()));
  }

  virtual void visit(FunctionSpaceMeshFunction& f)
  {
    //FIXME: implement me!
  }

  virtual void visit(FunctionSpaceUndefined& f)
  {
    CFD_EXCEPTION("Tried to evaluate undefined function space");
  }
};

}

}

#endif
