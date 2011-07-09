#ifndef SIMPLE_CFD_LAGRANGE_TRIANGLE_HPP
#define SIMPLE_CFD_LAGRANGE_TRIANGLE_HPP

#include <map>
#include <vector>
#include <utility>
#include <cassert>
#include <cmath>
#include <ostream>
#include <numeric>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include "simple_cfd_fwd.hpp"
#include "exception.hpp"
#include "mesh.hpp"
#include "finite_element.hpp"
#include "coordinate_transformation.hpp"
#include "nodal_basis_solver.hpp"
#include "numeric/tensor.hpp"
#include "numeric/tensor_size.hpp"
#include "numeric/index.hpp"
#include "numeric/excafe_expression.hpp"
#include "numeric/convert_expression.hpp"
#include "numeric/math_utilities.hpp"
#include "numeric/functional.hpp"
#include "cell_manager.hpp"
#include "exception.hpp"
#include "capture/assembly/position_placeholder.hpp"
#include "capture/assembly/scalar_placeholder.hpp"
#include "capture/assembly/generic_symbol.hpp"
#include "capture/assembly/tensor_operations.hpp"
#include "mp/rational.hpp"

namespace cfd
{

template<std::size_t R>
class LagrangeTriangle : public FiniteElement<2>
{
public:
  typedef TriangularCell cell_type;
  static const std::size_t rank = R;
  static const std::size_t dimension = cell_type::dimension;

  typedef Tensor<dimension> value_type;
  typedef Tensor<dimension> gradient_type;
  typedef Tensor<dimension> divergence_type;
  typedef vertex<dimension> vertex_type;

private:
  static const unsigned int tensor_size = detail::Power<dimension, rank>::value;
  typedef vertex<dimension, mp::Rational> exact_vertex_type;
  typedef boost::tuple<MeshEntity, exact_vertex_type, detail::ScalarPlaceholder::expression_t> dof_info_t;

  const cell_ref_t referenceCell;
  const std::size_t degree;
  std::vector<dof_info_t> dofInfo;

  std::vector< ExcafeExpression<detail::ScalarPlaceholder> > buildPrincipalBasis() const
  {
    using detail::ScalarPlaceholder;

    std::vector< ExcafeExpression<ScalarPlaceholder> > bases;

    detail::PositionPlaceholder position;
    const ScalarPlaceholder x = position[0];
    const ScalarPlaceholder y = position[1];
    const ScalarPlaceholder sf(detail::GenericSymbol("sf"));
    const ExcafeExpression<ScalarPlaceholder> sfExpr(sf);

    // These are the substitutions required to map coordinates from the
    // UFC reference triangle, to the the reference square.
    ExcafeExpression<ScalarPlaceholder>::value_map varChange;
    varChange.bind(x, 2*(2*x)/sf - 1);
    varChange.bind(y, 2*y-1);

    // We explicitly substitute in the scale factors separately so that
    // we can expand first, and avoid the problem of sfs in the
    // denominator not cancelling. 
    
    // The Sherwin-Karniadakis book differs from the 1994 paper on the
    // triangular spectral element method on the multiplier in the
    // scaling factor. We use the book, and have divided it by two.

    ExcafeExpression<ScalarPlaceholder>::value_map sfChange;
    sfChange.bind(sf, 1-y);

    for(std::size_t m=0; m<degree+1; ++m)
    {
      for(std::size_t n=0; m+n<degree+1; ++n)
      {
        const ExcafeExpression<ScalarPlaceholder> pm = 
          MathUtilities::jacobi(x, 0, 0, m).substituteValues(varChange);
        const ExcafeExpression<ScalarPlaceholder> pn = 
          MathUtilities::jacobi(y, 2*m+1, 0, n).substituteValues(varChange);

        ExcafeExpression<ScalarPlaceholder> basis = (pm * pn * pow(sfExpr, m)).normalised();
        basis = basis.substituteValues(sfChange).normalised();
        bases.push_back(basis);
      }
    }

    return bases;
  }

  std::vector<dof_info_t> buildDofInfo() const
  {
    const std::vector< ExcafeExpression<detail::ScalarPlaceholder> > principals = buildPrincipalBasis();
    const std::map< MeshEntity, std::vector<exact_vertex_type> > pointMap = referenceCell->getPoints(degree);

    typedef std::map< MeshEntity, std::vector<exact_vertex_type> >::value_type point_mapping_t;
    std::vector<exact_vertex_type> points;

    BOOST_FOREACH(const point_mapping_t& mapping, pointMap)
      points.insert(points.end(), mapping.second.begin(), mapping.second.end());

    NodalBasisSolver<dimension> nodalBasisSolver(points, principals);
    const std::vector< ExcafeExpression<detail::ScalarPlaceholder> > bases = nodalBasisSolver.getBases();

    std::vector<dof_info_t> dofInfo;
    std::size_t dof = 0;
    BOOST_FOREACH(const point_mapping_t& mapping, pointMap)
    {
      BOOST_FOREACH(const exact_vertex_type& v, mapping.second)
      {
        assert(dof < points.size());
        dofInfo.push_back(boost::make_tuple(mapping.first, v, bases[dof]));
        ++dof;
      }
    }

    return dofInfo;
  }

  static value_type evaluate(const tensor_expr_t& tensor, const vertex_type& v)
  {
    expression_t::value_map valueMap;
    detail::PositionPlaceholder position;

    for(std::size_t i=0; i<dimension; ++i)
      valueMap.bind(position[i], v[i]);

    ExpressionEvaluator<expression_t> evaluator(valueMap);
    value_type result(tensor.getSize());

    std::transform(tensor.begin(), tensor.end(), result.begin(), evaluator);
    return result;
  }

public:
  // We define the numbering of bases on a cell in the following fashion
  // index_into_tensor * number_of_nodes_on_cell + node_on_cell_id

  LagrangeTriangle(const std::size_t _degree) : 
  referenceCell(CellManager::getInstance<cell_type>()), degree(_degree), dofInfo(buildDofInfo())
  {
  }

  std::size_t getRank() const
  {
    return rank;
  }

  std::size_t getDimension() const
  {
    return dimension;
  }

  tensor_expr_t getBasis(const std::size_t i, const detail::PositionPlaceholder& v) const
  {
    assert(i < spaceDimension());

    const TensorSize tensorSize(rank, dimension);
    const TensorIndex tensorIndex = TensorIndex::unflatten(tensorSize, i / dofInfo.size(), row_major_tag());

    tensor_expr_t basis(tensorSize);
    basis[tensorIndex] = boost::get<2>(dofInfo[i % dofInfo.size()]);

    return basis;
  }

  value_type evaluateTensor(const std::size_t i, const vertex_type& vRef) const
  {
    detail::PositionPlaceholder position;
    const tensor_expr_t basis = getBasis(i, position);
    return evaluate(basis, vRef);
  }

  gradient_type evaluateGradient(const std::size_t i, const vertex_type& vRef) const
  {
    detail::PositionPlaceholder position;
    detail::TensorOperations<dimension> tensorOps;
    const tensor_expr_t gradientExpr = tensorOps.grad(getBasis(i, position));
    return evaluate(gradientExpr, vRef);
  }

  divergence_type evaluateDivergence(const std::size_t i, const vertex_type& vRef) const
  {
    detail::PositionPlaceholder position;
    detail::TensorOperations<dimension> tensorOps;
    const tensor_expr_t divergenceExpr = tensorOps.gradToDiv(tensorOps.grad(getBasis(i, position)));
    return evaluate(divergenceExpr, vRef);
  }

  unsigned spaceDimension() const
  {
    return dofInfo.size() * tensor_size;
  }

  vertex_type getDofCoordinateGlobal(const Mesh<dimension>& m, const cell_id cid, const std::size_t dof) const
  {
    const CellVertices<dimension> vertices = m.getCoordinates(cid);
    return referenceCell->referenceToPhysical(vertices, getDofCoordinateLocal(dof));
  }

  vertex_type getDofCoordinateLocal(const std::size_t i) const
  {
    const exact_vertex_type vExact =  boost::get<1>(dofInfo[i % dofInfo.size()]);
    vertex_type v;

    for(std::size_t i=0; i<dimension; ++i)
      v[i] = vExact[i].toDouble();

    return v;
  }

  // NOTE: by permitting mapping dofs to tensor indices, this commits
  // us to using standard bases.
  unsigned getTensorIndex(const Mesh<dimension>& mesh, const std::size_t cid, const std::size_t dof) const
  {
    assert(dof < spaceDimension());
    return dof / dofInfo.size();
  }

  virtual std::vector< std::set<dof_t> > resolveIdenticalDofs(const Mesh<dimension>& m, const MeshEntity& entity, const std::set<dof_t>& dofsOnEntity) const
  {
    typedef std::map<std::size_t, std::set<dof_t> > tensor_index_to_dofs_map;
    tensor_index_to_dofs_map tensorIndexToDofsMap;

    BOOST_FOREACH(const dof_t& dof, dofsOnEntity)
    {
      assert(dof.getElement() == this);
      tensorIndexToDofsMap[dof.getIndex() / dofInfo.size()].insert(dof);
    }

    std::vector< std::set<dof_t> > sharedDofs;
    BOOST_FOREACH(const tensor_index_to_dofs_map::value_type& indexMapping, tensorIndexToDofsMap)
      sharedDofs.push_back(indexMapping.second);

    return sharedDofs;
  }

  virtual std::set< Dof<dimension> > getDofsOnEntity(MeshTopology& topology, const cell_id cid, const MeshEntity& entity) const
  {
    const std::size_t localIndex = referenceCell->getLocalIndex(topology, cid, entity);
    const MeshEntity localEntity(entity.getDimension(), localIndex);

    std::set< Dof<dimension> > result;

    for(std::size_t i=0; i < spaceDimension(); ++i)
    {
      if (boost::get<0>(dofInfo[i % dofInfo.size()]) == localEntity)
        result.insert(Dof<dimension>(this, cid, i));
    }

    return result;
  }

  virtual cell_ref_t getCell() const
  {
    return referenceCell;
  }

  void write(std::ostream& o) const
  {
    o << "finite_element(name=\"Lagrange Triangle\", rank=" << rank << ", dimension=" << dimension;
    o << ", space=" << spaceDimension() << ")";
  }
};

}

#endif
